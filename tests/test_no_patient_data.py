"""No patient-level records may be tracked in this repository.

data/raw/harmonized_cohort/ once held 612 patient records, built in part
from AACR Project GENIE, GENIE BPC, MSK-IMPACT and MSK-CHORD, carrying sex,
age, race, ethnicity, smoking status, stage and survival. Those sources'
data use agreements do not permit redistribution, and the files were tracked
for publication while contributing nothing: the pipeline reproduces every
derived table byte for byte without them.

They are untracked and gitignored now. This test is what stops them, or
anything like them, from coming back in a later `git add -A`. It inspects
only the header line of each tracked table, never the rows.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

# Column names that mark a row as belonging to an identifiable individual.
# Deliberately narrow: variant-level tables legitimately carry things like
# "n_samples", and an over-broad rule would be switched off rather than fixed.
PATIENT_COLUMNS = {
    "patient_key", "patient_id", "patientid",
    "sample_id", "sampleid", "tumor_sample_barcode",
    "submitter_id", "case_id", "bcr_patient_barcode",
}


def _tracked_tables() -> list[Path]:
    """Every .csv and .tsv git is tracking. Untracked files cannot be pushed."""
    out = subprocess.run(
        ["git", "ls-files", "-z", "*.csv", "*.tsv"],
        cwd=ROOT, capture_output=True, text=True, check=True).stdout
    return [ROOT / p for p in out.split("\0") if p]


def test_repository_is_a_git_checkout():
    """Guard the guard: a silent empty file list would pass vacuously."""
    assert (ROOT / ".git").exists(), "not a git checkout, so this test proves nothing"
    assert _tracked_tables(), "no tracked tables found; git ls-files returned nothing"


@pytest.mark.parametrize("path", _tracked_tables(), ids=lambda p: str(p.name))
def test_tracked_table_has_no_patient_identifier(path: Path):
    if not path.exists():
        pytest.skip(f"{path.name} is tracked but absent from the working tree")
    with path.open(encoding="utf-8", errors="replace") as fh:
        header = fh.readline()
    sep = "\t" if path.suffix == ".tsv" else ","
    cols = {c.strip().strip('"').lower() for c in header.split(sep)}
    offending = sorted(cols & PATIENT_COLUMNS)
    rel = path.relative_to(ROOT)
    assert not offending, (
        f"{rel} is tracked by git and its header carries {offending}, which "
        "identifies rows as individual patients. Patient-level records must "
        "not be published: untrack the file and add it to .gitignore."
    )


# ---------------------------------------------------------------- metadata --
def test_deposit_metadata_agrees_with_itself():
    """.zenodo.json and CITATION.cff describe the same deposit.

    They are separate files read by different services, so they drift
    silently. A Zenodo record carrying a different author or licence from
    the repository it archives is worse than no record.
    """
    import json

    import yaml

    z = json.loads((ROOT / ".zenodo.json").read_text())
    c = yaml.safe_load((ROOT / "CITATION.cff").read_text())

    assert z["upload_type"] == "software"
    assert z["access_right"] == "open"
    assert z["license"] == c["license"], "licence differs between the two files"

    zc, cc = z["creators"], c["authors"]
    assert len(zc) == len(cc) == 1, "author lists differ in length"
    assert cc[0]["orcid"].endswith(zc[0]["orcid"]), "ORCID differs"
    assert zc[0]["name"] == f'{cc[0]["family-names"]}, {cc[0]["given-names"]}'
    assert zc[0]["affiliation"] == cc[0]["affiliation"], "affiliation differs"

    repo = "https://github.com/allkure/mutagen-reachability"
    assert c["repository-code"] == repo
    assert any(r["identifier"] == repo for r in z["related_identifiers"])
