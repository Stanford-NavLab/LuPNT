"""Helpers to host heavy interactive figures outside the repo and embed them in the docs.

The rendered documentation shows the stored outputs of each ``ex*.ipynb``
(nbsphinx_execute = "never"). Interactive Plotly 3-D plots and Cesium globes are large
(megabytes of embedded mesh/trace data), so instead of baking them into the notebooks we
write each to a standalone HTML file under ``build/doc_assets/`` and embed a small
``<iframe>`` pointing at a copy hosted on the ``lupnt-doc-assets`` GitHub Pages site.

Workflow:
    1. ``pixi run render-notebooks``  -> writes the HTML files into ``build/doc_assets/``
       and stores only the tiny iframe in each notebook.
    2. Publish ``build/doc_assets/`` to the ``lupnt-doc-assets`` repo (see
       ``scripts/publish_doc_assets.sh``) so the iframe URLs resolve.

Override the base URL with the ``LUPNT_DOC_ASSETS_URL`` environment variable if the assets
live elsewhere (e.g. a fork or a different host).
"""

from __future__ import annotations

import os
from pathlib import Path

from IPython.display import HTML

ASSETS_BASE = os.environ.get(
    "LUPNT_DOC_ASSETS_URL",
    "https://stanford-navlab.github.io/lupnt-doc-assets",
).rstrip("/")


def _repo_root() -> Path:
    for base in [Path.cwd(), *Path.cwd().parents]:
        if (base / "python" / "pylupnt" / "__init__.py").exists():
            return base
    return Path.cwd()


STAGE = _repo_root() / "build" / "doc_assets"


def _embed(url: str, height: int, label: str) -> HTML:
    return HTML(
        f'<p style="margin:0 0 6px">&#9654; <a href="{url}" target="_blank" '
        f'rel="noopener"><b>Open the interactive {label} &#8599;</b></a></p>'
        f'<iframe src="{url}" width="100%" height="{height}" frameborder="0" '
        f'loading="lazy" allowfullscreen></iframe>'
    )


def embed_plotly(fig, name: str, height: int = 640, label: str = "3-D view") -> HTML:
    """Write ``fig`` to ``build/doc_assets/plots/<name>.html`` (Plotly.js from CDN) and
    return a small iframe pointing at the hosted copy. Replaces ``fig.show()``."""
    out = STAGE / "plots"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{name}.html"
    fig.write_html(str(path), include_plotlyjs="cdn", full_html=True)
    url = f"{ASSETS_BASE}/plots/{name}.html"
    print(f"[doc-asset] plots/{name}.html  ({path.stat().st_size / 1024:.0f} kB)  ->  {url}")
    return _embed(url, height, label)


def embed_cesium_scene(scene, name: str, height: int = 520, label: str = "globe") -> HTML:
    """Save a ``CesiumScene`` to ``build/doc_assets/cesium/<name>.html`` and return a small
    iframe pointing at the hosted copy. Replaces ``scene.show(...)``."""
    out = STAGE / "cesium"
    out.mkdir(parents=True, exist_ok=True)
    path = scene.save(out / f"{name}.html")
    url = f"{ASSETS_BASE}/cesium/{name}.html"
    print(f"[doc-asset] cesium/{name}.html  ({path.stat().st_size / 1024:.0f} kB)  ->  {url}")
    return _embed(url, height, label)
