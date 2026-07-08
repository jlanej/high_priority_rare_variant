# Idea doc — the CRAM-access phase

**Status:** proposal / not scheduled. A dependency-bundling note, not a build spec.

Several high-value roadmap items are individually "low effort to add here" but share one **non-trivial
new data-model dependency**: re-accessing the Kids First source **CRAMs** (and, for some, a
≥100–150-sample batch). Today the pipeline is deliberately VCF-native — Steps 1–7 never open an
alignment; only **Step 8** (igv mini-CRAMs) touches CRAMs, and only to slice tiny review windows.
This doc groups the CRAM-gated work so the access cost is paid once, in a coherent phase, instead of
piecemeal.

See **[docs/ROADMAP.md](ROADMAP.md)** for the ranked gap list this draws from.

## Why bundle

Re-accessing controlled-access CRAMs is the expensive, slow, approval-gated step — not the analysis
on top of it. Each item below is cheap *given* aligned reads and pointless *without* them. Bundling
means one data-access request, one bind-mount contract, one batch, and one round of provenance/QC.

The bridge already exists: `resources.cram_map` (TSV `sample<TAB>cram_path`) is defined today for the
Step 8 export. The phase reuses that exact contract as its single input — no new config surface for
the common case.

## What CRAM access unlocks

| Capability | What the CRAMs give that the VCF can't | Roadmap | Notes |
|---|---|---|---|
| **somalier** relate + ancestry | Genotypes at ~17.5k common polymorphic sites (our rare-variant VCFs are far too sparse to fingerprint a sample) | #2 / spine #1 | Cross-cohort **duplicate/swap** detection (a proband under two IDs fabricates recurrence — the one failure mode nothing currently catches) + per-sample **ancestry PCs**. Needs somalier sites VCF + 1KG/HGDP `--labels`. |
| **verifyBamID2** true FREEMIX | Read-level contamination estimate (upgrades Step 0's VCF-only CHARR fallback to the primary path) | #3 (done, VCF-only) | Step 0 already ingests `*.selfSM` if supplied; this phase is where those files get **produced**. |
| **GATK-gCNV** germline CNV | Read-depth over intervals across a batch → copy-number calls + ACMG/ClinGen dosage (AnnotSV/ClassifyCNV) | #9 | Largest true blind spot (single-exon RB1/SMARCB1/DICER1/NF1/PMS2 deletions). **Requires the batch**, not just per-trio reads. Feed CNV-in-trans into the comp-het resolver. |
| **Read-backed phasing** (WhatsHap) | Direct phase when both variants of a candidate comp-het pair lie within one read/fragment | #11 | Complements trio phasing for the parent-of-origin-only cis pairs that currently inflate carrier counts. Pairs naturally with #9's CRAM access. |
| **Pre-refinement PL/GT cross-check** (top hits) | (Uses the pre-refinement VCF, not the CRAM, but IGV-level read review of the locus is the human backstop) | gotcha | The gnomAD-prior-suppression failure mode; CRAM windows make the manual cross-check reviewable in the Step 8 igv export. |

Adjacent but **out of this phase** (different data, not just aligned germline reads): the somatic
second-hit / LOH overlay needs matched **tumor** CRAMs (PBTA/OpenPedCan), a separate access request.

Not CRAM-gated (stays VCF-native, do independently): ROH (`bcftools roh`), UPD screen (UPDhmm),
PP1/BS4 co-segregation, UTRannotator, calibrated recurrence null (already done).

## Downstream unlock: CoCoRV

Ancestry PCs from somalier are the prerequisite for **CoCoRV** — ancestry-stratified external-control
(gnomAD) burden with an empirical-null λ and discrete-aware FDR. That is the natural phase-2 after the
calibrated recurrence null (#1, done): it moves the headline "recurrent + constrained" signal from an
internal binomial tail to a case-vs-control burden test. CoCoRV is *gated on this phase* and does not
itself need CRAMs — only the ancestry labels this phase produces.

## Sketch (ordering within the phase)

1. **Access + bind contract.** Provision CRAMs; populate `resources.cram_map`; verify readability +
   reference concordance (`samtools quickcheck`, `@SQ` md5s vs the pipeline GRCh38). Fail loudly on
   any missing/mismatched member — same discipline as trio resolution.
2. **Per-sample fingerprint QC (somalier).** `extract` → `relate` across **all** resolved samples
   (not per-family) for dup/swap/relatedness; `ancestry` for PCs. Fold results into the Step 0 gate
   and a new cohort-level QC artifact. This is the highest value-per-effort item — do it first.
3. **verifyBamID2** over the batch → `*.selfSM` → Step 0's existing FREEMIX path (retires the CHARR
   fallback for these samples).
4. **GATK-gCNV** batch calling + dosage annotation → CNV candidate track; wire CNV-in-trans into the
   comp-het resolver and the gene consolidation.
5. **Read-backed phasing** for comp-het pairs the trio can't phase.
6. **(Phase 2) CoCoRV** using the ancestry labels from step 2.

## Cost / requirements

- Controlled-access data-use approval for the CRAMs (public repo, so **nothing here changes the
  no-PHI / `${ENV}`-paths golden rule** — CRAM paths stay in the git-ignored `cram_map`).
- A **batch** (≥100–150 samples) for GATK-gCNV to model read-depth; somalier/verifyBamID work
  per-sample but benefit from cohort context for relatedness.
- Extra images/resources: somalier sites VCF + 1KG labels, GATK gCNV models, verifyBamID2 resource
  panel. Bind-mounted, never baked in (same rule as the VEP cache).

## Open questions

- Batch composition/size for gCNV given the trios are not jointly genotyped — per-family vs pooled
  read-depth cohorts.
- Whether somalier cross-cohort dedup should **hard-block** recurrence counting or just flag (leaning
  flag + review, consistent with the never-drop ethos).
- Reference-build/interval-list parity between the CRAMs and the pipeline GRCh38.
