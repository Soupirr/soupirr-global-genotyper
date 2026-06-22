"""Shared paths and palette."""

import os

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA_FOLDER = os.path.join(_ROOT, "data")
SEQ_FOLDER = os.path.join(DATA_FOLDER, "sequences")
LOCATION_FOLDER = os.path.join(DATA_FOLDER, "locations")
PATHO_FOLDER = os.path.join(DATA_FOLDER, "pathogenicity")
HOSTS_FOLDER = os.path.join(DATA_FOLDER, "hosts")
MISC_FOLDER = os.path.join(_ROOT, "misc")
TOOLS_FOLDER = os.path.join(_ROOT, "tools")
VALIDATION_FOLDER = os.path.join(DATA_FOLDER, "validation")

# palette importé depuis ColorBrewer (Spectral 11 + Set1/Dark2 extensions)
PALETTE = [
    "#9e0142",
    "#d53e4f",
    "#f46d43",
    "#fdae61",
    "#d4a017",
    "#a8b400",
    "#5fb300",
    "#3d9e6b",
    "#66c2a5",
    "#3288bd",
    "#5e4fa2",
    "#e41a1c",
    "#377eb8",
    "#4daf4a",
    "#984ea3",
    "#ff7f00",
    "#a65628",
    "#f781bf",
    "#b15928",
    "#cab2d6",
]

CUSTOM_CSS = """
    <style>
    h1, h2, h3 {
        letter-spacing: 3px;
        font-weight: 200;
        background: linear-gradient(90deg, #00c9a7, #0099cc, #6699cc);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        filter: drop-shadow(0 0 20px rgba(0, 201, 167, 0.2));
    }
    header { visibility: hidden; }
    section[data-testid="stMain"] {
        background: radial-gradient(ellipse at top, #0c1a24 0%, #060d14 70%);
    }

    /* Plotly toolbar transparent */
    .js-plotly-plot .plotly .modebar-container {
        background: rgba(0,0,0,0) !important;
    }
    .js-plotly-plot .plotly .modebar-group {
        background: rgba(0,0,0,0) !important;
    }
    .js-plotly-plot .plotly .modebar {
        background: rgba(0,0,0,0) !important;
    }
    .js-plotly-plot .plotly .modebar-btn {
        background: rgba(0,0,0,0) !important;
    }
    .js-plotly-plot .plotly .modebar-btn path {
        fill: rgba(0, 201, 167, 0.4) !important;
    }
    .js-plotly-plot .plotly .modebar-btn:hover path {
        fill: rgba(0, 201, 167, 0.9) !important;
    }
    .stProgress > div > div > div > div {
        background-color: #00aebc;  /* red */
    }
    [data-testid="stSelectboxVirtualDropdownEmpty"] {
        display: none !important;
    }
"""
