# Supplementary reachability track

`reachability_track.html` is the interactive per-residue reachability
inventory. It renders the same data as Figure 3C, one gene and one drug at a
time, at full residue resolution: pick a gene, pick an inhibitor, hover a
residue to read its top-ranked substitution, its forced substitution class,
and each mutagen's one-step supply for that class.

It is a lookup view, not a test. Nothing in it is a p-value.

## Files

| file | what it is |
|---|---|
| `reachability_track.html` | the page, one self-contained file |
| `reachability_track.json` | its data, a byte copy of `data/derived/reachability_track.json` |
| `reachability_<gene>_<drug>.svg` | one vector export per gene-drug pair, 7 in all |
| `reachability_<gene>_<drug>.pdf` | the same seven, as PDF |

The page is written to the Artifact contract, so it carries no `<html>`,
`<head>` or `<body>` of its own: the platform supplies the document skeleton
at publish time. The export script below wraps it so the same single file
also renders in a plain browser.

## The data file is a copy

`reachability_track.json` here is a copy of the pipeline's output. Copies go
stale, so after any `make numbers` that touches the track, re-copy it:

```bash
cp data/derived/reachability_track.json paper/supplement/
```

`tests/test_manuscript_tables.py::test_supplement_track_matches_derived`
fails if you forget. It was added after this copy sat at the old
1,205-variant KRAS universe for eleven days while the pipeline had moved to
981.

## Regenerating the SVG and PDF exports

GenomeSpy Core 1.0 exports real SVG through its `imageExport.svg()`, so the
seven figures are vector rather than screenshots. The export is scripted, not
done by hand, so it cannot drift from the data file.

The page loads GenomeSpy from jsDelivr when opened normally. The export runs
in a container whose network allowlist does not reach that host, so it serves
a local copy of the bundle instead and passes it to the page through the
`?gs=` override. The bundle is third-party and is not committed; fetch it:

```bash
# 1. the GenomeSpy bundle, into paper/supplement/gs/ (gitignored)
npm pack @genome-spy/core@1.0.0
tar xzf genome-spy-core-1.0.0.tgz
mkdir -p paper/supplement/gs && cp -r package/dist/bundle/* paper/supplement/gs/

# 2. Playwright, for the headless render
npm i playwright        # Chromium is preinstalled at /opt/pw-browsers
```

Then, from the repository root:

```bash
node scripts/export_reachability_svg.mjs
```

It prints one line per export with the SVG's size, its rectangle count, and
the number of rasterized runs. A correct run reports `0 rasters` and
`0 warnings` on every line: any raster run means part of the figure fell back
to a bitmap and the result is no longer a vector figure.

The script writes SVG only. For the PDFs:

```bash
cd paper/supplement
for f in reachability_*.svg; do cairosvg "$f" -o "${f%.svg}.pdf"; done
```

Neither step is part of `make all`, and nothing in the manuscript's number
register depends on either. Both need network access and a browser, which the
reproduction path deliberately does not.
