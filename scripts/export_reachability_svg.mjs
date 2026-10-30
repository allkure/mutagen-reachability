/**
 * Render paper/supplement/alk_reachability.html headless and export one
 * vector SVG per ALK drug.
 *
 * GenomeSpy Core 1.0 exports real SVG (its `imageExport.svg`), so the
 * supplementary figure is vector, not a screenshot. The export is scripted
 * rather than done by hand so the figure cannot drift from
 * data/derived/reachability_track.json.
 *
 * The page loads GenomeSpy from jsDelivr by default. This container's network
 * allowlist does not reach that host, so we serve a local copy of the bundle
 * and pass it through the page's `?gs=` development override.
 *
 *   node scripts/export_reachability_svg.mjs
 *
 * Requires: npm i playwright  (Chromium is preinstalled at /opt/pw-browsers)
 * and a local GenomeSpy bundle; see SUPPLEMENT.md for the two commands.
 */
import { chromium } from "playwright";
import http from "node:http";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const SITE = path.join(ROOT, "paper", "supplement");
const OUT = SITE;
const PORT = 8732;
// Filled from the data file: every gene and each of its drugs.

const TYPES = { ".html": "text/html", ".json": "application/json",
                ".js": "text/javascript", ".mjs": "text/javascript" };

// Minimal static server. The bundle lives under paper/supplement/gs/ and is
// not committed; SUPPLEMENT.md says how to put it there.
const server = http.createServer((req, res) => {
  const rel = decodeURIComponent(req.url.split("?")[0]).replace(/^\/+/, "");
  const file = path.join(SITE, rel || "reachability_track.html");
  if (!file.startsWith(SITE) || !fs.existsSync(file)) {
    res.writeHead(404).end("not found");
    return;
  }
  // alk_reachability.html is an artifact page: it carries no document
  // skeleton, because the Artifact platform supplies one at publish time.
  // Wrap it here so the same single source file renders in a plain browser.
  if (path.extname(file) === ".html") {
    res.writeHead(200, { "Content-Type": "text/html" });
    res.end('<!doctype html><html lang="en"><head><meta charset="utf-8">'
      + '<meta name="viewport" content="width=device-width,initial-scale=1">'
      + "</head><body>" + fs.readFileSync(file, "utf8") + "</body></html>");
    return;
  }
  res.writeHead(200, { "Content-Type": TYPES[path.extname(file)] ?? "application/octet-stream" });
  fs.createReadStream(file).pipe(res);
});
await new Promise((r) => server.listen(PORT, "127.0.0.1", r));

if (!fs.existsSync(path.join(SITE, "gs", "index.es.js"))) {
  console.error("missing paper/supplement/gs/index.es.js — see SUPPLEMENT.md");
  process.exit(1);
}

const browser = await chromium.launch({
  args: ["--enable-unsafe-swiftshader", "--use-gl=angle",
         "--use-angle=swiftshader", "--ignore-gpu-blocklist", "--no-sandbox"],
});
const page = await browser.newPage({ viewport: { width: 1280, height: 760 } });
const problems = [];
page.on("pageerror", (e) => problems.push(e.message));
page.on("console", (m) => m.type() === "error" && problems.push(m.text()));
// This container cannot reach fonts.googleapis.com, so the web-font stylesheet
// fails here and the page falls back to its system stack. GenomeSpy draws its
// own labels either way, so that failure does not affect the exported figure.
page.on("requestfailed", (r) => {
  if (!/fonts\.(googleapis|gstatic)\.com/.test(r.url())) {
    problems.push(`request failed: ${r.url()}`);
  }
});
const isEnvironmental = (s) => /ERR_TUNNEL_CONNECTION_FAILED|Failed to load resource/.test(s);

await page.goto(
  `http://127.0.0.1:${PORT}/reachability_track.html?gs=./gs/index.es.js`,
  { waitUntil: "networkidle" });
await page.waitForFunction(() => !!window.__api, null, { timeout: 30000 });

const plan = await page.evaluate(() => {
  const d = window.__data;
  return d.genes.map((g) => ({ gene: g, drugs: d.by_gene[g].drugs }));
});

const pick = async (host, label) => {
  const b = page.locator(`#${host} button`, { hasText: label });
  await b.first().click();
};

for (const { gene, drugs } of plan) {
  await pick("genes", gene);
  await page.waitForTimeout(2500);
  for (const drug of drugs) {
  await pick("drugs", drug);
  await page.waitForTimeout(1200);
  const r = await page.evaluate(async () => {
    const res = await window.__api.imageExport.svg({
      rasterization: { maxVectorInstances: 20000 },
    });
    return { text: await res.blob.text(), warnings: res.warnings,
             rasterized: res.rasterized.length };
  });
  const slug = `${gene}_${drug}`.toLowerCase().replace(/[^a-z0-9]+/g, "_");
  const out = path.join(OUT, `reachability_${slug}.svg`);
  fs.writeFileSync(out, r.text);
  const images = (r.text.match(/<image/g) ?? []).length;
  if (r.rasterized || images) {
    problems.push(`${slug}: ${r.rasterized} raster run(s), ${images} <image>`);
  }
  for (const w of r.warnings) problems.push(`${slug}: ${w}`);
  console.log(`${slug.padEnd(24)} ${(r.text.length / 1024).toFixed(0).padStart(4)} KB  `
    + `${(r.text.match(/<rect/g) ?? []).length} rects  `
    + `${images} rasters  ${r.warnings.length} warnings`);
  }
}

await browser.close();
server.close();

const real = problems.filter((p) => !isEnvironmental(p));
if (real.length) {
  console.error("\nexport is not clean:");
  for (const p of real) console.error("  " + p);
  process.exit(1);
}
console.log("\nevery export is pure vector, no warnings");
