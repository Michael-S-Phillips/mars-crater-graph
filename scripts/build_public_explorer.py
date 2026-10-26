"""Build the public, account-free versions of the explorer from explorer/page.html.

  python scripts/build_public_explorer.py --src explorer --docs docs \
      --standalone dist/mars_crater_graph_explorer.html \
      [--doi 10.5281/zenodo.XXXXXXX] [--repo https://github.com/<user>/mars-crater-graph]

Outputs
  explorer/index.html            full page with the tiled basemap (explorer/tiles/);
                                 opens from disk or any web server
  index.html (repository root)   redirect to explorer/, for GitHub Pages served from the root
  dist/...explorer.html          one self-contained file (data and coarse relief inlined;
                                 uses the tiles only when --tiles-url points at a hosted copy)

Tiles are cut by scripts/make_tiles.py where the source mosaics live. Levels missing
from a host (for example THEMIS level 5, kept out of the repository for size) are
detected by the page, which then falls back to the level below.

The full pages load their data with <script src="data/mcg_data.js">, which browsers
allow from file:// URLs; fetch() of JSON is blocked there, which is why the page
fragment alone fails when opened from disk.
"""
import argparse, base64, html, json, os, shutil

ap = argparse.ArgumentParser()
ap.add_argument("--src", required=True)
ap.add_argument("--docs", default=None, help="(unused; kept for older command lines)")
ap.add_argument("--standalone", required=True)
ap.add_argument("--doi", default="")
ap.add_argument("--repo", default="")
ap.add_argument("--tiles-url", default="", help="base URL of hosted tiles for the single-file build")
a = ap.parse_args()

body = open(os.path.join(a.src, "page.html"), encoding="utf-8").read()
# the hosted-artifact build carries its own <title> and a host-specific meta tag
title = body.split("<title>", 1)[1].split("</title>", 1)[0]
body = body.replace(f"<title>{title}</title>\n", "").replace(
    '<meta name="claude:viewport-fit-sides" content="cover">\n', "")

links = []
if a.repo:
    links.append(f'Code and data: <a href="{html.escape(a.repo)}">{html.escape(a.repo)}</a>.')
if a.doi:
    links.append(f'Cite: <a href="https://doi.org/{html.escape(a.doi)}">doi:{html.escape(a.doi)}</a>.')
if links:
    body = body.replace(
        'Host units: Tanaka et al. (2014).</p>',
        'Host units: Tanaka et al. (2014). ' + " ".join(links) + '</p>')

desc = ("Interactive map of 83,370 Mars craters linked by superposition: click a crater to "
        "trace the craters known to be older and younger than it.")
head = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>{title}</title>
<meta name="description" content="{html.escape(desc)}">
<meta property="og:title" content="{title}">
<meta property="og:description" content="{html.escape(desc)}">
<style>
/* page reset (the hosted artifact supplies these itself) */
:root {{ padding-top: env(safe-area-inset-top, 0px); padding-bottom: env(safe-area-inset-bottom, 0px); box-sizing: border-box; height: 100%; }}
body {{ margin: 0; }}
img {{ max-width: 100%; }}
a {{ color: inherit; }}
</style>
"""
tail = "\n</body>\n</html>\n"

d = os.path.join(a.src, "data")
# data as a script file, so the page works from disk as well as from a server
with open(os.path.join(d, "mcg_data.js"), "w", encoding="utf-8") as fh:
    fh.write("window.MCG_DATA = {"
             f"nodes: {open(os.path.join(d, 'nodes.json')).read()},"
             f"edges: {open(os.path.join(d, 'edges.json')).read()},"
             f"unlinked: {open(os.path.join(d, 'unlinked.json')).read()},"
             f"meta: {open(os.path.join(d, 'meta.json')).read()}"
             "};\n")
def tiles_config(base):
    return {"base": base, "default": "themis", "layers": [
        {"id": "themis", "name": "THEMIS day IR", "maxZ": 5,
         "credit": "Basemap: THEMIS daytime infrared global mosaic, 100 m/pixel "
                   "(Christensen et al., 2004; Edwards et al., 2011)."},
        {"id": "relief", "name": "HRSC-MOLA shaded relief", "maxZ": 4,
         "credit": "Basemap: HRSC-MOLA blended shaded relief, 200 m/pixel (Fergason et al., 2018)."},
    ]}


def tiles_script(base):
    return f"<script>window.MCG_TILES = {json.dumps(tiles_config(base))};</script>\n"


loader = '<script src="data/mcg_data.js"></script>\n'
page = head + "</head>\n<body>\n" + loader + tiles_script("tiles") + body + tail

# 1) full page next to its data and tiles, plus a root redirect for GitHub Pages
with open(os.path.join(a.src, "index.html"), "w", encoding="utf-8") as fh:
    fh.write(page)
root = os.path.dirname(os.path.abspath(a.src))
with open(os.path.join(root, "index.html"), "w", encoding="utf-8") as fh:
    fh.write('<!doctype html>\n<meta charset="utf-8">\n<title>Mars Crater Superposition Graph</title>\n'
             '<meta http-equiv="refresh" content="0; url=explorer/">\n'
             '<a href="explorer/">Open the Mars crater superposition graph explorer</a>\n')
open(os.path.join(root, ".nojekyll"), "w").close()

# 2) single self-contained file
relief = "data:image/jpeg;base64," + base64.b64encode(open(os.path.join(d, "relief.jpg"), "rb").read()).decode()
inline = ("<script>window.MCG_DATA = {"
          f"nodes: {open(os.path.join(d, 'nodes.json')).read()},"
          f"edges: {open(os.path.join(d, 'edges.json')).read()},"
          f"unlinked: {open(os.path.join(d, 'unlinked.json')).read()},"
          f"meta: {open(os.path.join(d, 'meta.json')).read()},"
          f"relief: {json.dumps(relief)}"
          "};</script>\n")
os.makedirs(os.path.dirname(a.standalone) or ".", exist_ok=True)
with open(a.standalone, "w", encoding="utf-8") as fh:
    fh.write(head + "</head>\n<body>\n" + inline
             + (tiles_script(a.tiles_url.rstrip("/")) if a.tiles_url else "") + body + tail)
print("page:", os.path.join(a.src, "index.html"))
print("root redirect:", os.path.join(root, "index.html"))
print("standalone:", a.standalone, f"{os.path.getsize(a.standalone)/1e6:.1f} MB")
