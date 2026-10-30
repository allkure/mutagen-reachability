#!/usr/bin/env python3
"""Per-residue reachability inventory for ALK, KRAS G12C and MET.

The Conclusions call the reachability inventory the auditable deliverable:
which resistance variants a mutagen-and-drug pair cannot produce. An inventory
is a thing you look things up in, not a summary statistic, so this script emits
it per residue, for every (drug, mutagen) pair, in the form a positional track
can read.

For each residue we take the substitution the chemistry-free atlas ranks
highest there, and ask what each mutagen's one-step supply for that
substitution is. A residue where the best available substitution falls in a
class the chemistry cannot make is a blind spot at that position.

Scope, which the caption has to carry: the inventory runs over class-forced
variants only, that is, substitutions a single nucleotide change can make
whose every route falls in one substitution class. Within that set a mutagen
is blocked only when its spectrum gives that class zero, which is why ENU and
Mutazyme II come out at zero blocked residues. It does not mean they have no
blind spot. The larger one, variants needing two or more nucleotide changes,
is outside this set by construction and is reported separately in the
manuscript, numerically for KRAS and MET. The Oh ALK atlas contains no
two-step variants at all by prime-edit design, so the three genes are not
comparable on that axis and this view does not attempt it.

This is a lookup view, not a test. We make no claim that reachability clusters
along the gene; nothing here is a p-value.

Supply spectra are imported from the same source the benchmark uses, so the
track cannot drift from the numbers in the manuscript.

Writes data/derived/reachability_track.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.spectrum import load_spectra  # noqa: E402

CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
ALK_DECOMP = ROOT / "data" / "derived" / "oh_decomposition_all_drugs.csv"
KRAS_DECOMP = ROOT / "data" / "derived" / "feng_decomposition.csv"
MET_DMS = ROOT / "data" / "raw" / "met_estevam2024" / "met_dms_wide.csv"
MET_CODONS = ROOT / "data" / "reference" / "codons_met.yaml"
OUT = ROOT / "data" / "derived" / "reachability_track.json"

# Expected class-forced universe per gene, as a guard against a silent
# change in the classifier. ALK is verified against the Oh atlas construct
# (331 codons, all reproducing the atlas wild-type residue). KRAS was 1,205
# until the classifier was corrected: that figure counted variants such as
# R68L and V8M as one-step when the genetic code puts them two changes from
# the wild-type codon. MET has no published figure, so it is unconstrained.
EXPECTED_N = {"ALK": 1427, "KRAS G12C": 981}
TRANSCRIPT = {"ALK": "NM_004304.5", "KRAS G12C": "NM_033360.4",
              "MET": "NM_000245.2"}
# MET sigma is the ROSACE effect under a drug minus that study's DMSO arm,
# on the exon-14-deleted background: the same quantity, and the same
# background, that the manuscript's MET section uses.
MET_BACKGROUND = "exon14del"
MET_DRUGS = {"NVP-BVU972": "delta_NVP", "AMG 458": "delta_A458"}

# Same four spectra the top-K benchmark uses (recompute_topk_with_ci.py).
# ENU and Mutazyme II come from the reference file; the two synthetic
# chemistries are definitional, not measured.
SPECTRA = load_spectra(ROOT / "data" / "reference" / "spectra.yaml")
MUTAGENS: dict[str, dict[str, float]] = {
    "ENU": dict(SPECTRA["ENU_Arnold2012_nominal"].values),
    "Mutazyme II": dict(SPECTRA["MutazymeII_DoiSupp2022"].values),
    "CBE+ABE": {"C>A": 0.0, "C>G": 0.0, "C>T": 0.5,
                "T>A": 0.0, "T>C": 0.5, "T>G": 0.0},
    "Twist/PE": {c: 1.0 / 6.0 for c in CLASSES},
}
# Order the lanes from chemistry-free to most restricted.
LANE_ORDER = ["Twist/PE", "Mutazyme II", "ENU", "CBE+ABE"]


def load_met_variants() -> pd.DataFrame:
    """MET class-forced variants with their chemistry-free selection score.

    The other two genes arrive pre-classified in a decomposition table. MET
    has none, so the classification is done here with met_audit's own
    `classify`, against met_audit's own codon table, so the two agree by
    construction rather than by coincidence.
    """
    sys.path.insert(0, str(ROOT / "scripts"))
    from met_audit import classify  # noqa: E402

    codons = {int(k): v for k, v in
              yaml.safe_load(MET_CODONS.read_text())["codons"].items()}
    dms = pd.read_csv(MET_DMS)
    dms = dms[dms.background == MET_BACKGROUND]

    rows = []
    for drug, col in MET_DRUGS.items():
        d = dms.dropna(subset=[col])
        for r in d.itertuples():
            codon = codons.get(int(r.residue))
            if codon is None:
                continue
            kind, cls, _ = classify(codon, r.mutation)
            if kind != "class-forced":
                continue
            rows.append({"aa_change": r.aa_change, "residue": int(r.residue),
                         "forced_class": cls, "drug": drug,
                         "sigma": float(getattr(r, col))})
    var = pd.DataFrame(rows)

    # Both drug arms must classify the same variants the same way; the codon
    # table does not depend on which inhibitor was assayed.
    per_drug = var.groupby("drug").size()
    assert per_drug.nunique() == 1 or per_drug.std() / per_drug.mean() < 0.25, \
        f"MET arms disagree in size: {dict(per_drug)}"
    return var


def load_variants(gene: str) -> pd.DataFrame:
    """Class-forced variants with their chemistry-free selection score.

    sigma is a property of (variant, drug) and does not depend on the
    mutagen, so the several mutagen rows per variant carry the same value and
    collapse to one. The two genes' decomposition tables differ in shape: the
    ALK one is already class-forced and carries a residue column, while the
    KRAS one holds every variant and names the residue `position`.
    """
    if gene == "MET":
        return load_met_variants()
    if gene == "ALK":
        d = pd.read_csv(ALK_DECOMP)
    else:
        d = pd.read_csv(KRAS_DECOMP)
        d = d[d["class_forced"] == True].rename(columns={"position": "residue"})

    key = ["aa_change", "residue", "forced_class", "drug"]
    chk = d.groupby(key)["drug_selec"].nunique()
    assert chk.max() == 1, f"{gene}: drug_selec varies within a (variant, drug) group"

    var = (d[key + ["drug_selec"]]
           .drop_duplicates()
           .rename(columns={"drug_selec": "sigma"})
           .reset_index(drop=True))
    assert set(var.forced_class) <= set(CLASSES), sorted(set(var.forced_class))

    per_drug = len(var) // var.drug.nunique()
    if gene in EXPECTED_N:
        assert per_drug == EXPECTED_N[gene], (
            f"{gene}: expected {EXPECTED_N[gene]:,} class-forced variants, "
            f"got {per_drug:,}")
    return var


def build(gene: str) -> dict:
    var = load_variants(gene)
    per_drug = len(var) // var.drug.nunique()

    # The top-sigma substitution at each residue, per drug.
    top = (var.sort_values("sigma", ascending=False)
              .groupby(["drug", "residue"], as_index=False)
              .first())

    records = []
    for lane, spec in MUTAGENS.items():
        # The published spectra are quoted to three decimals and do not sum to
        # exactly 1 (ENU 1.001, Mutazyme II 1.021). They are used as given, as
        # ratios, so the tolerance here is loose by design.
        assert abs(sum(spec.values()) - 1.0) < 0.05, (lane, sum(spec.values()))
        assert set(spec) == set(CLASSES), (lane, sorted(spec))
        t = top.copy()
        t["mutagen"] = lane
        t["S"] = t.forced_class.map(spec).astype(float)
        t["reachable"] = t.S > 0
        records.append(t)
    trk = pd.concat(records, ignore_index=True)

    # CBE+ABE must forbid exactly four of the six classes.
    be = MUTAGENS["CBE+ABE"]
    assert sum(v == 0 for v in be.values()) == 4, be

    residues = sorted(top.residue.unique())
    out = {
        "gene": gene,
        "transcript": TRANSCRIPT[gene],
        "residue_min": int(min(residues)),
        "residue_max": int(max(residues)),
        "n_residues": len(residues),
        "n_class_forced_variants": int(per_drug),
        "drugs": sorted(top.drug.unique()),
        "lane_order": LANE_ORDER,
        "classes": CLASSES,
        "spectra": {k: {c: float(v[c]) for c in CLASSES}
                    for k, v in MUTAGENS.items()},
        "sigma": [
            {"residue": int(r.residue), "drug": r.drug,
             "sigma": round(float(r.sigma), 4),
             "aa_change": r.aa_change, "forced_class": r.forced_class}
            for r in top.itertuples()
        ],
        "track": [
            {"residue": int(r.residue), "drug": r.drug, "mutagen": r.mutagen,
             "aa_change": r.aa_change, "forced_class": r.forced_class,
             "sigma": round(float(r.sigma), 4), "S": round(float(r.S), 4),
             "reachable": bool(r.reachable)}
            for r in trk.itertuples()
        ],
    }

    print(f"{gene}  {out['transcript']}  residues {out['residue_min']}"
          f"--{out['residue_max']}  (n = {out['n_residues']})  "
          f"{out['n_class_forced_variants']:,} class-forced variants per drug")
    for lane in LANE_ORDER:
        sub = trk[trk.mutagen == lane]
        for drug in out["drugs"]:
            s = sub[sub.drug == drug]
            n = int((~s.reachable).sum())
            print(f"    {lane:12s} {drug:11s} {n:4d} / {len(s):4d} residues "
                  f"out of reach ({100*n/len(s):.0f}%)")
    print()
    return out


def main() -> None:
    genes = {g: build(g) for g in ("ALK", "KRAS G12C", "MET")}
    out = {
        "genes": list(genes),
        "lane_order": LANE_ORDER,
        "classes": CLASSES,
        "spectra": {k: {c: float(v[c]) for c in CLASSES}
                    for k, v in MUTAGENS.items()},
        "by_gene": genes,
    }
    OUT.write_text(json.dumps(out, separators=(",", ":")))
    rows = sum(len(g["track"]) for g in genes.values())
    print(f"wrote {OUT}  ({OUT.stat().st_size/1024:.0f} KB, {rows:,} track rows)")


if __name__ == "__main__":
    main()
