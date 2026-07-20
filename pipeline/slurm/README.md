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
