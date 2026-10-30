# Archiving a release on Zenodo

The manuscript's Availability statement cites a Zenodo DOI. This is how to
mint it, and what has already been prepared.

## Prepared in this repository

| file | what it does |
|---|---|
| `.zenodo.json` | the deposit's metadata. Zenodo reads it at release time, so the record gets the right title, author, ORCID, affiliation, licence and keywords instead of guessing from the repo name. |
| `CITATION.cff` | GitHub renders a "Cite this repository" button from it, and reference managers read it. |

Both name the same author, ORCID and licence, and a test checks they agree
so they cannot drift apart.

## What gets archived, and what does not

Zenodo archives the repository as GitHub serves it at the tagged commit, so
the deposit contains exactly the tracked files and nothing that `.gitignore`
excludes.

Deliberately excluded:

- **Patient-level clinical records** (`data/raw/harmonized_cohort/*.csv`).
  Built in part from AACR Project GENIE, GENIE BPC, MSK-IMPACT and
  MSK-CHORD, whose data use agreements do not permit redistribution, and
  carrying demographics and survival for 612 individuals. No result depends
  on them: the pipeline reproduces every derived table byte for byte
  without them, and `tests/test_no_patient_data.py` fails if any tracked
  table regains a patient identifier column.
- **The vendored GenomeSpy bundle** (`paper/supplement/gs/`), third-party
  code fetched from npm. See `paper/supplement/SUPPLEMENT.md`.

A Zenodo DOI is permanent and a record cannot be quietly withdrawn, so
check `git ls-files` before tagging if anything about the data has changed.

## Steps

1. **Push the repository** and confirm it is public. A private repository
   is invisible to Zenodo, and returns the same 404 to a reviewer as one
   that does not exist.

   ```bash
   git push -u origin main
   ```

2. **Link the repository.** Sign in to [zenodo.org](https://zenodo.org) with
   GitHub, open the GitHub settings page, find `allkure/mutagen-reachability`
   and switch it on. Do this *before* creating the release: Zenodo only
   archives releases made after the switch is flipped.

3. **Tag and release.**

   ```bash
   git tag -a v1.0.0 -m "Version archived with the manuscript"
   git push origin v1.0.0
   ```

   Then publish a GitHub release from that tag. Zenodo picks it up within a
   few minutes and mints two DOIs: one for this version, and a *concept DOI*
   that always resolves to the latest version.

4. **Cite the version DOI, not the concept DOI**, in the manuscript. A
   reader following the paper's link should land on the code as it was when
   the results were produced, not on a later revision.

5. **Fill the three placeholders** that are waiting on this:

   - `paper/paper.tex`, Availability of data and materials: replace
     `DOI\_TO\_BE\_INSERTED` with the version DOI.
   - `CITATION.cff`: uncomment and fill `doi` and `date-released`.
   - `README.md`, Citation section: link the preprint and the release.

   Then rebuild: `make paper`.

## Order of operations

Post the bioRxiv preprint first, then archive, then submit. The preprint
timestamps the work, and the Zenodo record can carry the preprint DOI as a
related identifier. If you archive first you will want to edit the record's
metadata afterwards to add it, which is possible but avoidable.
