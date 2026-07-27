# Distributed run on SLURM

Runs the pipeline as a coherent SLURM job graph, with **Step 2 (VEP) scattered one
contig per node** and the rest chained by dependencies:

```
prep ─afterok→ plan ─submits→ scatter[array 0..N] ─afterok→ gather ─afterok→ downstream
(Steps 0–1)    (enumerate       (one contig of        (concat shards +     (Steps 3–8)
               contigs, submit   VEP+CADD per task)    split-vep + output)
               the rest)
```

**Coherence.** Every edge is `--dependency=afterok`, so any failure halts everything
downstream — you never get a call set built from a partial Step 2. As a second guard,
`gather` independently re-verifies that every contig shard is complete and dies loudly
if one is missing. So neither the scheduler nor the pipeline alone can silently produce
a truncated result.

**Why this shape.** The scatter array's size is the number of contigs with variants,
which is only known after `prep` builds the union — so `plan` (a tiny job that runs
after `prep`) enumerates them and submits the array + gather + downstream. That is the
standard SLURM pattern for a dynamically-sized array; it needs no `sbatch --array=N`
guess and no over-provisioning.

## Use

```bash
cp cluster.env.example cluster.env
$EDITOR cluster.env            # partition/account, container launch, per-phase resources
./submit_slurm.sh cluster.env # run once, from a login node
```

> **Prerequisite: `HPRV_CONFIG` must be a fully-resolved config (literal paths, not `${ENV}`
> placeholders).** Every phase runs `apptainer exec --cleanenv`, which strips host env vars, so
> `config.py` inside the container cannot expand `${REF_FASTA}` / `${VEP_CACHE}` / `${CADD_*}` and
> `run_pipeline.sh` dies at preflight. Resolve the config once before submitting — e.g. source your
> resource env (`prepare_resources.sh … emit-env`) and `envsubst < config.example.yaml > config.yaml`.
> (An interactive `run_pipeline.sh` on a login node does *not* hit this, because it inherits your shell's env.)

`submit_slurm.sh` prints the `prep` and `plan` job IDs. The scatter/gather/downstream
IDs are chosen by `plan` at runtime — find them after `plan` runs in
`$HPRV_WORK/slurm_jobids.txt`, or with `squeue -u $USER --name=hprv-scatter,hprv-gather,hprv-down`.

## Resume / retry

The whole graph is idempotent through the pipeline's `.done` sentinels:

- **Walltime kill or a transient failure:** just re-run `./submit_slurm.sh cluster.env`.
  Completed pieces are skipped; only unfinished work re-runs. A failed contig means
  `gather` is left with `DependencyNeverSatisfied` — requeue the failed array element
  (`scontrol requeue <arrayjobid>_<idx>`) or re-submit; its shard has no `.done`, so
  only that contig re-runs, and `gather` proceeds once it completes.
- **Force a full re-annotation** (e.g. you changed the cache or CADD): remove the shard
  directory first — `rm -rf $HPRV_WORK/annotate_shards/` — then re-submit.

## Ending the graph early (e.g. run Step 8 elsewhere)

Step 8 (mini-CRAM / igv export) reads the **source CRAM store**, which at some sites is a
FUSE/network mount visible only on login/interactive nodes — not the batch compute nodes the
graph runs on. To stop the distributed graph at Step 7 and run Step 8 where the CRAMs live, set
in `cluster.env`:

```sh
DOWN_TO=7
```

The `downstream` job then runs Steps 3–7, and `submit_slurm.sh` finishes with the annotated call
set, gene ranking and xlsx in `$HPRV_WORK`. Run Step 8 afterwards on a node that can see the
CRAMs — it only reads files already in `$HPRV_WORK`, so it is a plain one-off:

```sh
apptainer exec --cleanenv --bind "$HPRV_BINDS" "$HPRV_SIF" \
    run_pipeline.sh --config "$HPRV_CONFIG" --from 8 --to 8
```

(`DOWN_FROM` is configurable too, for symmetry; the range must satisfy `3 <= FROM <= TO <= 8`.)
The audit summary is re-assembled at the end of each run, so it ends up reflecting Step 8 once
that runs.

## Sizing notes

- **`SCATTER_TIME`** is per *contig*, not per genome. Size it to the largest contig
  (chr1 ≈ 8% of a WGS union — a few hours), never to the whole 57M.
- **`SCATTER_CONCURRENCY`** (the `%K` in `--array=0-N%K`) is how many nodes Step 2
  occupies at once.
- **CADD staging.** The per-shard bottleneck is CADD's per-variant `tabix` reads into
  the ~81 GB file. If it lives on a network/FUSE mount, stage it to node-local disk at
  the top of the `scatter` phase (an `rsync`/`cp` prologue, then point `CADD_SNV` at the
  local copy) — 81 GB × concurrent-nodes of transfer only pays off if per-contig VEP
  runtime ≫ copy time, so measure one contig first.

## Correctness

The distributed path produces **byte-identical output to a single in-process VEP run** —
only the `vep` call is scattered; `split-vep` and every guard run once on the reassembled
whole. This is enforced by `tests/integration/assert_shard_equivalence.sh` (sharded ==
single) and exercised in CI.

## Not yet distributed

Steps 0/1/4 (per-trio, currently serial) are the next candidates for trio-arrays — on
WGS, Step 1 (the union build) becomes the tall pole once Step 2 is scattered. The `.done`
idempotency they already carry makes that the natural next increment; see
[docs/limitations.md](../../docs/limitations.md).

---

## Optional: Step 8b (NHF) as a per-trio array

Step 8b is the dominant cost after Step 2 on a large cohort — **classification, not kraken2 DB
load, dominates** (~20 min/trio on WGS), so 221 trios is ~74 h serial. It parallelises per trio
almost linearly.

```bash
# AFTER a completed Step 8 (the array reads the mini-CRAMs + per-trio VCFs that Step 8 produces)
sbatch phase.sbatch nhf-plan
```

`nhf-plan` writes `$WORK/igv/nhf_trios.txt` (only trios with outstanding work, so a resubmit
schedules a smaller array), submits `nhf-scatter` as a job array (one task = one trio), then a
dependent `nhf-gather` that folds the results into `variants.tsv`.

Knobs live in `cluster.env`: `NHF_CPUS` / `NHF_MEM` / `NHF_TIME`, `NHF_CONCURRENCY`, and
`NHF_GATHER_*`. Set `outputs.igv.nonhuman_screen.threads` to match `NHF_CPUS` — it is **decoupled**
from `outputs.igv.extract_jobs`, which bounds FUSE-gentle CRAM slicing and is deliberately small.

Three things that differ from the Step-2 scatter:

- **It needs Step 8 to have run.** Sub-tasks deliberately skip Pass 1/2, which rebuild state shared
  across trios (the per-trio VCF copy that 8b's content key hashes, the candidate BEDs, the extract
  task list). Concurrent rebuilds would mean torn reads and unstable keys.
- **No source-CRAM mount is touched.** 8b reads only `$WORK/igv/{crams,vcfs}/`, so it is safe on
  batch compute even at sites where Step 8's slicing must run on an interactive node
  (`DOWN_TO=7` + Step 8 elsewhere).
- **Gather depends `afterany`, not `afterok`, and never dies.** An unscreened member is legitimate
  here — `members: carriers` skips hom-ref parents by design, and a failed classify degrades to
  blank NHF columns (distinct from a real `0.0`). Gather reports counts and writes `variants.tsv`
  regardless; re-run `nhf-plan` to see what is genuinely outstanding.

Keep `NHF_CONCURRENCY` modest: every task does random reads over one shared kraken2 DB (hundreds of
GB, far too large to page-cache), so a wide array is real shared-storage load. Do **not** stage the
DB to node-local scratch — it exceeds typical local disk and the copy costs more than it saves.
