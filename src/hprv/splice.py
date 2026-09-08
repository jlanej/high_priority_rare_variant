"""SpliceAI event decomposition and `SpliceAI=` parsing (pure; no VCF I/O, so it is unit-testable).

SpliceAI reports FOUR delta scores — acceptor gain / acceptor loss / donor gain / donor loss — and,
for each, the offset in bp (relative to the variant) of the splice site it predicts is created or
destroyed. The screen (`selection.py`) and Step 9's tier read only the MAX of the four. That is the
right quantity for a gate, and the wrong one for a reviewer, who needs to know WHICH event, WHERE,
and what it implies for the transcript. This module turns the components into that, once, so Step 5
(precomputed scores from the per-trio VCF INFO) and Step 5b (wide-window live rescoring) describe an
event in exactly the same vocabulary.

The ``effect`` vocabulary:

  site_loss                 one loss event above the floor: an annotated site is lost. Exon skipping
                            or intron retention; the FRAME needs the exon length, which this module
                            does not have (that is the GTF-backed follow-up).
  site_gain                 one gain event: a cryptic site, or one end of a pseudoexon. Frame unknown.
  cryptic_shift_in_frame    a gain AND a loss of the SAME site type (two donors, or two acceptors):
  cryptic_shift_frameshift  the exon boundary moves by |dp_gain - dp_loss| nt — in frame iff that is
                            a multiple of 3. This is computable with no exon structure at all, and it
                            is the most common shape of a cryptic-site activation.
  cryptic_shift             the same pair but a position is missing, so the frame cannot be computed.
  paired_gains              an acceptor gain AND a donor gain: a pseudoexon candidate.
  paired_losses             an acceptor loss AND a donor loss: a whole-exon-loss candidate.
  complex                   a gain and a loss of DIFFERENT site types.

The floor for calling an event is the screen's own ``filters.functional.spliceai_ds_min`` (0.2): the
same number that decides a variant HAS a splice signal decides which components count as events. A
variant whose max delta is below the floor gets no event at all (its max still rides along as
``spliceai_ds``), because naming the largest of four near-zero numbers an "event" would be noise
dressed as a finding.
"""

from __future__ import annotations

from typing import Optional

EVENTS = ("acceptor_gain", "acceptor_loss", "donor_gain", "donor_loss")
# SpliceAI's own suffixes, in the order the `SpliceAI=` string and the VEP plugin's CSQ keys use.
EVENT_SUFFIX = {"acceptor_gain": "AG", "acceptor_loss": "AL", "donor_gain": "DG", "donor_loss": "DL"}
SITE = {e: e.split("_")[0] for e in EVENTS}     # acceptor | donor
KIND = {e: e.split("_")[1] for e in EVENTS}     # gain | loss

EMPTY = {"event": None, "ds": None, "dp": None, "pos": None,
         "event2": None, "event2_ds": None, "event2_dp": None, "event2_pos": None,
         "shift_nt": None, "shift_frame": None, "effect": None}


def _num(x) -> Optional[float]:
    if x is None:
        return None
    s = str(x).strip()
    if s in ("", ".", "NA", "nan", "None"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _int(x) -> Optional[int]:
    v = _num(x)
    return None if v is None else int(v)


def decompose(ds, dp, pos=None, floor: float = 0.2) -> dict:
    """Name the event(s) behind four SpliceAI delta scores.

    ``ds`` / ``dp`` map event -> delta score / bp offset (missing values allowed). ``pos`` is the
    variant's 1-based position, so ``pos`` in the result is the genomic coordinate of the affected
    site (the place to look in RNA and in IGV). Ties in the delta score break in EVENTS order, so
    the result is deterministic.
    """
    ds = {e: _num((ds or {}).get(e)) for e in EVENTS}
    dp = {e: _int((dp or {}).get(e)) for e in EVENTS}
    out = dict(EMPTY)
    ranked = sorted((e for e in EVENTS if ds[e] is not None),
                    key=lambda e: (-ds[e], EVENTS.index(e)))
    if not ranked or ds[ranked[0]] < floor:
        return out

    def site_pos(e):
        return (pos + dp[e]) if (pos is not None and dp[e] is not None) else None

    e1 = ranked[0]
    out.update(event=e1, ds=ds[e1], dp=dp[e1], pos=site_pos(e1))
    e2 = next((e for e in ranked[1:] if ds[e] >= floor), None)
    if e2 is None:
        out["effect"] = "site_loss" if KIND[e1] == "loss" else "site_gain"
        return out
    out.update(event2=e2, event2_ds=ds[e2], event2_dp=dp[e2], event2_pos=site_pos(e2))
    kinds, sites = {KIND[e1], KIND[e2]}, {SITE[e1], SITE[e2]}
    if len(sites) == 1 and kinds == {"gain", "loss"}:
        gain, loss = (e1, e2) if KIND[e1] == "gain" else (e2, e1)
        if dp[gain] is None or dp[loss] is None:
            out["effect"] = "cryptic_shift"
            return out
        # The new site replaces the lost one, so the exon boundary moves by this many nt and the
        # exon length changes by the same amount (either strand: the MAGNITUDE is what sets frame).
        shift = dp[gain] - dp[loss]
        out["shift_nt"] = shift
        out["shift_frame"] = "in_frame" if shift % 3 == 0 else "frameshift"
        out["effect"] = "cryptic_shift_in_frame" if shift % 3 == 0 else "cryptic_shift_frameshift"
    elif kinds == {"gain"}:
        out["effect"] = "paired_gains"
    elif kinds == {"loss"}:
        out["effect"] = "paired_losses"
    else:
        out["effect"] = "complex"
    return out


def parse_spliceai_info(value) -> Optional[dict]:
    """Parse the live tool's ``SpliceAI=`` INFO value into ONE coherent entry.

    The format is ``ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL``, comma-joined
    across the genes overlapping the variant. The entry with the LARGEST max delta score is kept
    whole (first wins on a tie) so its scores and offsets stay together: a per-event max across
    genes — what the Step-2b backfill does, correctly, for a plain keep-path gate — would pair a
    score from one gene with an offset from another and mis-place the event.

    Returns ``{"symbol", "ds": {event: float|None}, "dp": {event: int|None}, "n_genes"}`` or None
    when nothing parsed.
    """
    if not value:
        return None
    best, n = None, 0
    for entry in str(value).split(","):
        parts = entry.split("|")
        if len(parts) < 6:
            continue
        n += 1
        ds = {e: _num(parts[2 + i]) for i, e in enumerate(EVENTS)}
        dp = {e: (_int(parts[6 + i]) if len(parts) > 6 + i else None) for i, e in enumerate(EVENTS)}
        present = [v for v in ds.values() if v is not None]
        if not present:
            continue
        cand = {"symbol": parts[1].strip(), "ds": ds, "dp": dp, "_max": max(present)}
        if best is None or cand["_max"] > best["_max"]:
            best = cand
    if best is None:
        return None
    best.pop("_max")
    best["n_genes"] = n
    return best


def max_ds(ds) -> Optional[float]:
    """The standard SpliceAI delta score: the max over the four components, or None."""
    vals = [_num(v) for v in (ds or {}).values()]
    vals = [v for v in vals if v is not None]
    return max(vals) if vals else None
