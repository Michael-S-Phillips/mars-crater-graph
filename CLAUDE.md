# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`mcgraph` builds a directed superposition (relative-age) graph of Mars impact craters from the Robbins and Hynek (2012) catalog (2020 release), plus a static HTML/JS map explorer. `README.md` is the authoritative description of the method, results, and validation. Its numbers (edge counts, agreement percentages) are taken from `outputs/summary.json` and `outputs/validation.json`; if the pipeline is rerun, update the README, `CITATION.cff`, and `.zenodo.json` to match.

## Python environment

The base Anaconda Python (`~/anaconda3/bin/python3`) cannot import scipy (missing `libgfortran.5.dylib`). The `sentinel` conda env runs the tests but has no `pyarrow`, so pipeline steps that read or write parquet need an env with everything in `requirements.txt`. `rasterio` is optional (`GeoTiffDEM` and the chip-extraction scripts only).

## Commands

Tests are plain functions run as a script (pytest is not required):

```bash
~/anaconda3/envs/sentinel/bin/python tests/test_core.py
# single test
~/anaconda3/envs/sentinel/bin/python -c "import sys; sys.path.insert(0,'tests'); import test_core as t; t.test_frac_circle_inside()"
```

The full pipeline is in the README under "Running the pipeline". It runs in order: `run_evidence.py` → `build_graph.py` → `validate.py` → `make_figures.py` → `export_explorer.py` → (`make_tiles.py`) → `build_public_explorer.py`. The pipeline needs external inputs that are not in the repo: the Robbins catalog CSV, the MOLA MEGDR 128 ppd tiles, the SIM 3292 geology shapefile, and the THEMIS and HRSC-MOLA mosaics for tiles. Scripts import `mcgraph` by inserting the repo root into `sys.path`; there is no installable package.

## Architecture and data flow

- `work/` (gitignored) holds intermediates: `catalog.parquet`, `pairs_evidence.parquet`, and the `mola128_global.i2` memmap mosaic built once from the 16 MEGDR tiles.
- `mcgraph/pairs.py` finds neighbors within 3 r of the larger crater using a cKDTree on unit vectors, then classifies each pair as `contain` / `overlap` / `touch` / `ejecta_zone`. Pairs are stored as `big`/`small` node indices.
- `mcgraph/evidence.py` orients `overlap` pairs. **Vote sign convention: +1 = the smaller crater is younger, −1 = the larger crater is younger, 0 = undecided.** E1 = MOLA rim-prominence continuity (`topo_*`), E2 = `ARC_IMG` rim-arc completeness (`arc_*`). E3 = `DEG_RIM` degradation is used only for validation and must never orient edges. `combine()` merges the votes and sets the confidence.
- `mcgraph/graph.py`: **edges point older → younger.** `edge_table(pairs, min_conf)` filters by confidence (`high` > `medium` > `low`). The default is `medium`, which excludes topography-only edges (`low`), possible basin rings, and near-concentric duplicates. It also computes `break_cycles`, node metrics (`depth`, `height`, transitive counts, components), and the Hasse diagram (transitive reduction).
- Node indices are row positions in the catalog. Public outputs are also keyed by Robbins `CRATER_ID`; keep both columns when changing exports.
- `mcgraph/units.py`: host geologic units are attributes only, never edges (see README "Limitations").

### Explorer build

- `explorer/page.html` is the **source** page. `scripts/export_explorer.py` writes `explorer/data/{nodes,edges,unlinked,meta}.json` and `relief.jpg` from `outputs/`.
- `scripts/build_public_explorer.py` (standard library only) generates `explorer/index.html`, `explorer/data/mcg_data.js`, the root `index.html` redirect, `.nojekyll`, and `dist/mars_crater_graph_explorer.html`, which has the data inlined. Edit `page.html` and rebuild; do not hand-edit the generated files.
- The built pages load data through `<script src="data/mcg_data.js">` because `fetch()` of JSON is blocked from `file://`. That is why `page.html` alone does not run from disk.
- Tiles: `explorer/tiles/{themis,relief}/<z>/...jpg`, cut by `make_tiles.py`. THEMIS level 5 (~750 MB) is local only. The page detects missing levels and falls back to the level below.

## Release and publishing constraints

`CLAUDE_CODE_RELEASE_PROMPT.md` is the step-by-step plan for publishing to GitHub, GitHub Pages, and Zenodo; follow it, including its **[ask]** checkpoints, when doing release work. Constraints that apply to any git work here:

- The folder sits on an external drive (`/Volumes/Fangorn`), where macOS writes `._*` and `.DS_Store` files. Never commit them, and set `git config core.fileMode false` in every repo.
- `main` holds the code, outputs, and built explorer **without** `explorer/tiles/`. The website (page, `mcg_data.js`, `relief.jpg`, and tile levels 0–4) goes on a separate `gh-pages` branch, built in a working copy on the internal disk and pushed in batches of ≤ ~250 MB.
- Do not use Git LFS (Pages does not serve LFS files). Check that no staged file exceeds 50 MB (GitHub rejects files over 100 MB).
- `explorer.zip` (~1.4 GB) is a hand-made archive and not part of the release.
- Do not delete, move, or rename files unless a release step says to. Never force-push or rewrite remote history.
- Once the DOI and repository URL exist, rebuild with `build_public_explorer.py --repo ... --doi ... --tiles-url <site>/explorer/tiles` so the footers and the single-file build link to them.
