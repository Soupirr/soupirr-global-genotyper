"""Shared paths, palette and CSS."""

import os

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA_FOLDER = os.path.join(_ROOT, "data")
SEQ_FOLDER = os.path.join(DATA_FOLDER, "sequences")
LOCATION_FOLDER = os.path.join(DATA_FOLDER, "locations")
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
    "#377eb8",
    "#4daf4a",
    "#984ea3",
    "#ff7f00",
    "#a65628",
    "#f781bf",
    "#b15928",
    "#cab2d6",
]

# Palette de couleur sans conotation sémantique pour les stats
PATHO_PALETTE = [
    "#A855F7",
    "#0BF5F1",
    "#0609D4",
    "#EC4899",
    "#FF006F",
]

# Look de l'appli
_CSS_COMMON = """
    [data-testid="stAppDeployButton"] { display: none !important; }
    [data-testid="stSelectboxVirtualDropdownEmpty"] { display: none !important; }
    [data-testid="stStatusWidget"] { display: none !important; }


    /* Plotly toolbar transparent */
    .js-plotly-plot .plotly .modebar-container,
    .js-plotly-plot .plotly .modebar-group,
    .js-plotly-plot .plotly .modebar,
    .js-plotly-plot .plotly .modebar-btn {
        background: rgba(0,0,0,0) !important;
    }
    .stProgress > div > div > div > div {
        background-color: #00aebc;
    }
"""

_CSS_DARK = """
    h1 {
        letter-spacing: 3px;
        font-weight: 600 !important;
        background: linear-gradient(90deg, #00c9a7, #0099cc, #6699cc);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        filter: drop-shadow(0 0 20px rgba(0, 201, 167, 0.2));
    }
    h2, h3 {
        letter-spacing: 3px;
        font-weight: 200;
        background: linear-gradient(90deg, #00c9a7, #0099cc, #6699cc);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        filter: drop-shadow(0 0 20px rgba(0, 201, 167, 0.2));
    }
    section[data-testid="stMain"] {
        background: radial-gradient(ellipse at top, #0c1a24 0%, #060d14 70%);
    }
    .js-plotly-plot .plotly .modebar-btn path {
        fill: rgba(0, 201, 167, 0.4) !important;
    }
    .js-plotly-plot .plotly .modebar-btn:hover path {
        fill: rgba(0, 201, 167, 0.9) !important;
    }
"""

_CSS_LIGHT = """
    h1 {
        letter-spacing: 3px;
        font-weight: 600 !important;
        background: linear-gradient(90deg, #007a66, #0066aa, #3355aa);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        filter: drop-shadow(0 0 10px rgba(0, 122, 102, 0.15));
    }
    h2, h3 {
        letter-spacing: 3px;
        font-weight: 200;
        background: linear-gradient(90deg, #007a66, #0066aa, #3355aa);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        filter: drop-shadow(0 0 10px rgba(0, 122, 102, 0.15));
    }
    .js-plotly-plot .plotly .modebar-btn path {
        fill: rgba(0, 122, 102, 0.5) !important;
    }
    .js-plotly-plot .plotly .modebar-btn:hover path {
        fill: rgba(0, 122, 102, 1.0) !important;
    }
    /* Plotly text (axes, legends, titles) en sombre sur fond clair */
    .js-plotly-plot .plotly text {
        fill: #1a2e3d !important;
    }
"""


def get_custom_css(theme: str = "Dark") -> str:
    theme_css = _CSS_LIGHT if theme == "Light" else _CSS_DARK
    return f"<style>{_CSS_COMMON}{theme_css}</style>"


CUSTOM_CSS = get_custom_css("Dark")
