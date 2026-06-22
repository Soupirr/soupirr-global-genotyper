"""
Newcastle Disease Virus F Gene Genotyper
Streamlit Web Application
"""

import streamlit as st
from ndv_genotyper.tabs import map_tab, tree_tab, analyze_tab, stats_tab, validation_tab
from ndv_genotyper.config import CUSTOM_CSS

# Page configuration
st.set_page_config(
    page_title="NDV Genotyper",
    page_icon="misc/icon.png",
    layout="wide",
    initial_sidebar_state="collapsed",
)
st.title("Newcastle Disease Virus F Gene Genotyper")
st.markdown("""
Analyze Newcastle Disease Virus sequences to identify genotypes and pathogenicity.
""")
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)

# read the hidden flag from the URL (?dev=1)
show_validation = st.query_params.get("dev") == "1"

# def des onglets
tab_labels = [
    "Analyze Sequences",
    "Phylogenetic Trees",
    "Global Distribution Map",
    "Help - Statistics",
]
if show_validation:
    tab_labels.append("Precision Validation")

tabs = st.tabs(tab_labels)
tab_analyze, tab_tree, tab_map, help_tab = tabs[0], tabs[1], tabs[2], tabs[3]

# ============================================================================
# TAB 1: ANALYZE SEQUENCE
with tab_analyze:
    analyze_tab.render()
references = analyze_tab.load_all_references()[0]


# TAB 2: PHYLOGENETIC TREES
with tab_tree:
    tree_tab.render(references)


# TAB 3: GLOBAL DISTRIBUTION MAP
with tab_map:
    map_tab.render()


# TAB 4: Stats/Help
with help_tab:
    stats_tab.render()
# ============================================================================

# DEV TAB: VALIDATION
if show_validation:
    with tabs[4]:
        validation_tab.render(references)
