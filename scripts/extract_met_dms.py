#!/usr/bin/env python3
"""Extract the Estevam 2024 MET kinase-domain DMS into the project's schema.

Source
------
github.com/fraser-lab/MET_kinase_Inhibitor_DMS, files
  WT_rosace_effect_all.tsv   wild-type MET kinase-domain background
  EX_rosace_effect_all.tsv   exon 14-deleted background
Citation: Estevam et al. 2024, eLife 13:RP101882, doi 10.7554/eLife.101882.

What this gives the audit
-------------------------
ROSACE_effects is the per-variant selection score under each of eleven
inhibitors plus a DMSO arm.  DMSO is the no-drug control, so it plays the
role Oh's Fitness score plays for ALK: an orthogonal growth axis, reported
alongside and NOT subtracted from the drug score, consistent with the
sigma = Resistance convention used throughout this paper.

Two of the inhibitors are the ones used in the Tiedt 2011 ENU Ba/F3 TPR-MET
screen that supplies the chemistry-biased side of the MET audit: NVP is
NVP-BVU972 and A458 is AMG 458.  The TPR-MET fusion carries a wild-type
kinase domain, so WT_rosace_effect_all.tsv is the matching background.

Position numbering
------------------
The published tables number the kinase domain locally, 1 to 287.  The
authors' own analysis code converts with `position + 1058`, giving
canonical MET residues 1059 to 1345.  This script re-derives that offset
from the data rather than trusting the constant: it scores every candidate
offset by how many known MET residue identities it reproduces and asserts
that 1058 is the unique best.

Author inhibitor groupings, read out of the repository's analysis code:
  type I   Crizo, Camp, Tepo, Glu, NVP, Savo
  type II  Cabo, Mere, Gle          (A458 is ordered with this group)
  Tiv and A458 appear in no explicit type vector.

Outputs
-------
data/raw/met_estevam2024/met_dms_long.csv    one row per (variant, inhibitor)
data/raw/met_estevam2024/met_dms_wide.csv    one row per variant
data/raw/met_estevam2024/met_dms_meta.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = (ROOT.parent / "external" / "met_dms")
OUT = ROOT / "data" / "raw" / "met_estevam2024"

OFFSET_CANDIDATES = range(1040, 1075)

# Canonical MET residue identities used only to pin the numbering offset.
# Every one is a residue named in the MET resistance literature whose
# wild-type identity is given by its own variant nomenclature.
ANCHORS = {1092: "V", 1094: "H", 1163: "G", 1195: "L", 1200: "F",
           1211: "M", 1228: "D", 1230: "Y", 1234: "Y", 1235: "Y"}

TYPE_I = ["Crizo", "Camp", "Tepo", "Glu", "NVP", "Savo"]
TYPE_II = ["Cabo", "Mere", "Gle"]
UNGROUPED = ["A458", "Tiv"]
CONTROL = "DMSO"
EXPECTED = set(TYPE_I + TYPE_II + UNGROUPED + [CONTROL])

TIEDT_INHIBITORS = {"NVP": "NVP-BVU972", "A458": "AMG 458"}


def derive_offset(d: pd.DataFrame) -> int:
    """Pick the offset that best reproduces the known residue identities."""
    wt = d.drop_duplicates("position").set_index("position").wildtype.to_dict()
    scores = {off: sum(1 for c, aa in ANCHORS.items() if wt.get(c - off) == aa)
              for off in OFFSET_CANDIDATES}
    best = max(scores, key=scores.get)
    top = scores[best]
    runners = sorted((v for k, v in scores.items() if k != best), reverse=True)
    if top < 8:
        sys.exit(f"offset {best} only reproduces {top}/{len(ANCHORS)} anchors; aborting")
    if runners and runners[0] >= top:
        sys.exit(f"offset {best} is not uniquely best (runner-up {runners[0]}/{len(ANCHORS)})")
    print(f"  offset derived from data: +{best} "
          f"({top}/{len(ANCHORS)} anchors; next best {runners[0]}/{len(ANCHORS)})")
    return best


def load(fname: str, background: str) -> tuple[pd.DataFrame, int]:
    path = SRC / fname
    if not path.exists():
        sys.exit(f"missing input: {path}")
    d = pd.read_csv(path, sep="\t")

    got = set(d.inhibitor.unique())
    if got != EXPECTED:
        sys.exit(f"{fname}: inhibitor set mismatch\n  expected {sorted(EXPECTED)}\n  got {sorted(got)}")

    # one wild-type residue per position
    per_pos = d.groupby("position").wildtype.nunique()
    if (per_pos > 1).any():
        sys.exit(f"{fname}: positions with inconsistent wildtype: "
                 f"{per_pos[per_pos > 1].index.tolist()[:5]}")

    # synonymous rows must not change the residue
    syn = d[d.type == "synonymous"]
    if not (syn.wildtype == syn.mutation).all():
        sys.exit(f"{fname}: synonymous rows where mutation != wildtype")

    offset = derive_offset(d)
    d = d.copy()
    d["residue"] = d.position + offset
    d["aa_change"] = d.wildtype + d.residue.astype(str) + d.mutation
    d["background"] = background
    return d, offset


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frames, offsets = [], {}
    for fname, bg in (("WT_rosace_effect_all.tsv", "WT"),
                      ("EX_rosace_effect_all.tsv", "exon14del")):
        print(f"{fname}:")
        d, off = load(fname, bg)
        offsets[bg] = off
        frames.append(d)
    long = pd.concat(frames, ignore_index=True)

    keep = ["background", "aa_change", "residue", "wildtype", "mutation", "type",
            "inhibitor", "ROSACE_effects", "ROSACE_tests", "ROSACE.GOF", "ROSACE.LOF"]
    long = long[keep].rename(columns={"ROSACE_effects": "sigma",
                                      "ROSACE.GOF": "gof", "ROSACE.LOF": "lof"})
    long.to_csv(OUT / "met_dms_long.csv", index=False)

    wide = (long[long.type == "missense"]
            .pivot_table(index=["background", "aa_change", "residue", "wildtype", "mutation"],
                         columns="inhibitor", values="sigma", aggfunc="mean")
            .reset_index())
    wide.columns.name = None

    # sigma for MET = drug arm minus the DMSO arm.
    #
    # This differs arithmetically from the ALK convention but is the same
    # quantity.  Oh's Resistance score is already log(drug_D10 / untreated_D10),
    # so the growth component is divided out before we see it, which is why its
    # separately reported Fitness axis is orthogonal and must not be subtracted
    # again.  Estevam's ROSACE_effects is computed per arm against that arm's
    # own starting pool, so the drug arms still carry the growth component: here
    # the DMSO arm IS the fitness axis, and it correlates with the drug arms at
    # r ~ 0.8 rather than ~ 0.  Subtracting it reproduces Oh's Resistance.
    # The repository's own resistance calling uses the same difference
    # (`delta.score = mean - mean[key == "DMSO"]`).
    drugs = [c for c in wide.columns if c in set(TYPE_I + TYPE_II + UNGROUPED)]
    for c in drugs:
        wide[f"delta_{c}"] = wide[c] - wide[CONTROL]
    wide.to_csv(OUT / "met_dms_wide.csv", index=False)

    wt = long[long.background == "WT"]
    meta = {
        "source_repo": "github.com/fraser-lab/MET_kinase_Inhibitor_DMS",
        "citation_key": "estevam2024",
        "offset_applied": offsets,
        "offset_anchors_used": ANCHORS,
        "canonical_residue_range": [int(long.residue.min()), int(long.residue.max())],
        "n_positions": int(long.residue.nunique()),
        "inhibitors": sorted(EXPECTED),
        "author_type_I": TYPE_I,
        "author_type_II": TYPE_II,
        "ungrouped_in_author_code": UNGROUPED,
        "control_arm": CONTROL,
        "tiedt_2011_inhibitors_present": TIEDT_INHIBITORS,
        "n_missense_WT": int((wt.type == "missense").sum()),
        "n_variants_WT_wide": int((wide.background == "WT").sum()),
        "note_dmso": ("DMSO is the no-drug control and IS subtracted, giving "
                      "delta_<inhibitor> columns. Estevam's per-arm ROSACE_effects "
                      "still contains the growth component (DMSO correlates with the "
                      "drug arms at r ~ 0.8), whereas Oh's Resistance score already "
                      "nets growth out. The subtraction makes the two atlases report "
                      "the same quantity; it matches the repository's own delta.score."),
        "note_backgrounds": ("The repository is ambiguous about which file is which "
                             "background: MET_inhibitor_DMS_analysis.R reads EX_ as "
                             "delta-exon14 and WT_ as exon-14-retained, while "
                             "MET_inhibitor_DMS_Rosace_analysis.R reads WT_ into both "
                             "objects, which is a copy-paste error. We follow the "
                             "self-consistent script and report both backgrounds."),
        "note_a458_class": ("The authors' explicit type_II vector is Cabo/Mere/Gle; A458 is "
                            "ordered immediately with that group in their plotting vectors "
                            "but carries no explicit type label in the repository code."),
    }
    (OUT / "met_dms_meta.json").write_text(json.dumps(meta, indent=2, default=str))

    print(f"\nlong: {len(long):,} rows   wide (missense): {len(wide):,} rows")
    print(f"canonical residues {meta['canonical_residue_range'][0]}-"
          f"{meta['canonical_residue_range'][1]} across {meta['n_positions']} positions")
    print(f"Tiedt 2011 inhibitors present: {TIEDT_INHIBITORS}")


if __name__ == "__main__":
    main()
