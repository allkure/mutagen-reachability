# Chemical reachability limits mutagen-based TKI-resistance screens

Analysis code and derived data for the manuscript of the same name.

**Michele Ronco** — European Commission, Joint Research Centre, Ispra
[ORCID 0000-0002-2160-2452](https://orcid.org/0000-0002-2160-2452)

---

## What this is about

Preclinical Ba/F3 mutagenesis screens are a first line of evidence for
predicting kinase-inhibitor resistance, and every mutagen used in them has a
non-uniform substitution spectrum. This work asks whether that chemistry
limits *which* resistance variants a screen can find, or merely reorders the
ones it does.

It separates a screen's expected output into a supply term and a selection
term,

```
P(v)  ∝  S_mutagen(v) · exp(β·σ(v))
```

where `S` comes from codon tables and the mutagen's published six-class
spectrum, and `σ` comes from a chemistry-free deep mutational scan. The two
halves are then tested separately against three kinases — KRAS G12C, ALK and
MET — seven inhibitors and four mutagen chemistries.

The headline result is that chemistry acts as a **reachability threshold, not
a ranking weight**. A cytosine-and-adenine base-editor pair produces only two
of the six substitution classes, and across all three kinases it cannot reach
the best available substitution at 63%–66% of residues. ENU and Mutazyme II,
which supply all six classes at a non-zero rate, cannot reach none.

---

## Quick start

```bash
git clone https://github.com/allkure/mutagen-reachability.git
cd mutagen-reachability
pip install -r requirements.txt
make all          # derived tables, figures, both audits, tests
```

Or open the notebook, which redraws every figure and prints the manuscript's
numbers beside the register that produced them:

```bash
jupyter notebook notebooks/reproduce_figures.ipynb
```

Individual targets:

| command | what it does |
|---|---|
| `make numbers` | recompute the derived tables and the number register |
| `make figures` | redraw Figures 1–5 from the derived tables |
| `make audit` | check the register against the manuscript, both directions |
| `make test` | run the test suite |
| `make paper` | build `paper/paper.pdf` (needs `pdflatex` and `bibtex`) |

---

## Every number comes from a script

This is the part worth knowing about before reading anything else.

Each quantitative claim in the manuscript is computed by a named function and
recorded in **`data/derived/paper_numbers.json`**, rendered for humans as
**`paper/NUMBERS.md`**. Each entry carries the value, the formula, the inputs,
the function that produced it, and where it appears in the paper.

Two scripts check that register against the manuscript, in opposite
directions:

- **`scripts/audit_paper_numbers.py`** — every registered value must appear in
  `paper.tex`, in the section the register claims. It splits the file at
  `\section{Background}`, so a value registered for the Abstract that sits
  only in Results is caught.
- **`scripts/audit_unregistered_numbers.py`** — every number written in the
  prose must be traceable to a registered value, a documented constant, or a
  table cell the test suite verifies against its source table.

The second one is the one that earns its keep. It found six stale recovery
figures and four stale variant ranks shortly before submission, all left
behind by a codon-table rebuild that moved the ALK class-forced universe from
1,259 to 1,427 variants: the registered numbers moved with it, and the ones
living only in prose did not.

`tests/test_manuscript_tables.py` extends the same idea to the hand-typed
LaTeX tables, checking all 96 cells of Table 2 against `benchmark_r3.csv` and
Table 4 against `sigma_transport.json`. Rows are matched on their labels, not
their order, because the manuscript deliberately omits two rows the CSV
contains. A later correction to the KRAS codon classifier moved six of those
cells, and the test named each one.

The same pass turned up the limit of a register that stores only ranges. The
six per-anchor permutation *p*-values were quoted individually in the prose
but registered only as a summary, so one of them sat at a stale `0.78` while
the range it belonged to still matched. They are registered individually now.

Both audits and the full test suite run in CI on every push.

On a fixed environment the pipeline is deterministic: running `make all`
twice reproduces every table under `data/derived/` byte for byte, and
`data/derived/CHECKSUMS.sha256` is unchanged between runs. Only the two
generated Markdown files carry a timestamp in their header, so those are the
only files a repeat run rewrites.

---

## Layout

```
├── paper/                  manuscript, figures and the number register
│   ├── paper.tex           the manuscript
│   ├── NUMBERS.md          every registered number, human-readable
│   ├── PIPELINES.md        the three-kinase audit grid
│   ├── figure_[1-5].py     one script per figure, shared style in _style.py
│   └── supplement/         interactive reachability track and vector exports
├── scripts/                the pipeline, in dependency order in the Makefile
├── src/                    importable modules: codon tools, spectra, stats
├── tests/                  97 tests, including the table and data checks
├── notebooks/              reproduce_figures.ipynb, generated not hand-edited
├── data/
│   ├── raw/                as extracted from the primary papers, immutable
│   ├── reference/          codon tables and mutagen spectra
│   └── derived/            everything the pipeline computes
├── preregistration.md      Phase 0, locked 2026-10-13, before any analysis
├── ZENODO.md               how the archived release is minted
├── CITATION.cff            machine-readable citation
└── .zenodo.json            deposit metadata, read by Zenodo at release time
```

The notebook is generated by `scripts/build_reproduction_notebook.py`. Edit
that script, not the `.ipynb`.

---

## Data sources

The derived tables are redistributed with back-references to their sources.
The primary data remain under the terms of their original publications; see
[`data/LICENSE`](data/LICENSE).

| what | source |
|---|---|
| KRAS G12C saturation mutagenesis atlas | Feng et al. 2022 |
| ALK prime-editing resistance atlas | Oh et al. 2026 |
| MET kinase-domain deep mutational scan | Estevam et al. 2024 |
| MET ENU screen hit lists and clone counts | Tiedt et al. 2011, Tables S1–S2 |
| Germline ENU six-class spectrum | Arnold et al. 2012 |
| Mutazyme II spectrum | Doi et al. 2022, supplementary NGS |

`data/derived/CHECKSUMS.sha256` fixes the version of every derived table used
for the manuscript.

### What is deliberately not here

Patient-level clinical records are not distributed. They are built in part
from AACR Project GENIE, GENIE BPC, MSK-IMPACT and MSK-CHORD, whose data use
agreements do not permit redistribution, and they carry demographics and
survival for 612 individuals. No result in the manuscript depends on them:
with all three files removed, `make all` completes, both audits pass, and
every derived table reproduces byte for byte.
`data/raw/harmonized_cohort/README.md` keeps their schema, and
`tests/test_no_patient_data.py` fails if any tracked table regains a patient
identifier column.

---

## Environment

Python 3.13.16, with versions pinned in `requirements.txt` and recorded by
`scripts/env_versions.py` into `data/derived/env_versions.json`, which §2.14
of the paper quotes. The pinned set is the one the figures and numbers were
produced with.

`scripts/export_reachability_svg.mjs` is the one exception: it drives a
headless browser to export the supplementary reachability track as vector SVG,
and needs Node, `playwright`, and a local copy of the GenomeSpy bundle.
[`paper/supplement/SUPPLEMENT.md`](paper/supplement/SUPPLEMENT.md) gives the
commands. Nothing in `make all` depends on it.

---

## Licence

Code is MIT ([`LICENSE`](LICENSE)). Derived data are CC BY 4.0
([`data/LICENSE`](data/LICENSE)). Primary data keep the terms of their
original publications.

---

## Citation

A preprint and an archived release will be linked here once they are posted.
[`ZENODO.md`](ZENODO.md) describes how the release is archived and which
placeholders it fills.
