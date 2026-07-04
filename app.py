"""Streamlit Web Application - GENOTYPER"""

import streamlit as st
from genotyper.tabs import tree_tab, analyze_tab, stats_tab, validation_tab
from genotyper.config import CUSTOM_CSS, SEQ_FOLDER, MISC_FOLDER
from genotyper.migration import migrate_fasta_text
import json
import os
import csv

# Page configuration
st.set_page_config(
    page_title="Soupirr's Genotyper",
    page_icon="misc/icon.png",
    layout="wide",
    initial_sidebar_state="collapsed",
)
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)

# def des onglets
tab_labels = [
    "Analyze Sequences",
    "Phylogenetic Trees",
    "Statistics",
    "Precision Validation",
]


# Charge les config de pathogénécité et les motifs depuis les dossier d'entrée
def load_entry_config(entry_path: str) -> dict:
    config_path = os.path.join(entry_path, "_config.json")
    if not os.path.exists(config_path):
        return {}  # retourne rien si il n'y en a pas
    config = {}

    with open(config_path) as f:
        raw = json.load(f)
    if "cleavage_start" in raw:
        config["cleavage_start"] = raw["cleavage_start"]

    entry_name = os.path.basename(entry_path)
    motifs_path = os.path.join(entry_path, f"{entry_name}_motifs.csv")
    if os.path.exists(motifs_path):
        motifs_by_type = {}
        with open(motifs_path, newline="") as f:
            for row in csv.DictReader(f):
                motif = row["motif"].strip().upper()
                label = row["label"].strip()
                type_name = row["type"].strip().lower()
                if type_name not in motifs_by_type:
                    motifs_by_type[type_name] = {}
                motifs_by_type[type_name][motif] = label
        if motifs_by_type:
            config["motifs_by_type"] = motifs_by_type

    return config


st.html("""
    <style>

        /* Cache le menu kebab (trois petits points) en haut à droite */
        #MainMenu, [data-testid="stMainMenu"] {
            visibility: hidden;
        }

        /* Top bar transparente pour se fondre avec le fond de l'app */
        [data-testid="stHeader"] {
            background-color: transparent;
        }
    </style>
""")

# ============================================================================
# SIDEBAR

if "form_key" not in st.session_state:
    st.session_state["form_key"] = 0


def display_documentation():
    with open(os.path.join(MISC_FOLDER, "QUICK_START.md"), "r", encoding="utf-8") as f:
        st.markdown(f.read())


def reset_app():
    st.cache_resource.clear()
    st.cache_data.clear()
    st.session_state.pop("data_loaded", None)
    st.write("Cache cleared !")


with st.sidebar:
    documentation_button = st.button("Documentation", on_click=display_documentation)
    entry = sorted(os.listdir(SEQ_FOLDER))
    st.divider()
    selection = st.selectbox(
        "Select an entry", entry, index=None, placeholder="Select an entry..."
    )
    st.write("")
    with st.expander("Add new references datasets"):
        entry_name = st.text_input(
            "Enter the name of the entry", key=f"name_{st.session_state['form_key']}"
        )
        upload_file = st.file_uploader(
            "Upload the reference FASTA file",
            accept_multiple_files=True,
            type=["fasta", "fas", "fa", "txt"],
            key=f"fasta_{st.session_state['form_key']}",
        )

        st.divider()
        st.markdown(
            "**Pathogenicity configuration** *(optional - leave blank to skip)*"
        )

        cleavage_start = st.number_input(
            "bp before the virulence motif area",
            min_value=0,
            step=1,
            help="0-indexed nucleotide position where the cleavage/virulence site starts in the sequence.",
            key=f"cleavage_start_{st.session_state['form_key']}",
        )

        motif_file = st.file_uploader(
            "Upload virulence motif file",
            type=["csv", "txt"],
            help="CSV with columns: motif, label, type - where type is 'virulent' or 'avirulent'",
            key=f"csv_{st.session_state['form_key']}",
        )
        st.caption(
            "Expected format: `motif,label,type` - e.g. `RRQKRF,VFcs-1,virulent`"
        )

        adding_button = st.button("Add to the registry")
        if adding_button:
            if not upload_file:
                st.error("Please select at least one FASTA file")
            elif not entry_name.strip():
                st.error("Please enter a name")
            else:
                entry_path = os.path.join(SEQ_FOLDER, entry_name.strip())
                os.makedirs(entry_path, exist_ok=True)

                mig_reports = []
                for uploaded_files in upload_file:
                    dest_path = os.path.join(entry_path, uploaded_files.name)
                    raw_text = (
                        uploaded_files.getbuffer()
                        .tobytes()
                        .decode("utf-8-sig", errors="ignore")
                    )
                    migrated_text, mig_stats = migrate_fasta_text(raw_text)
                    with open(dest_path, "w", encoding="utf-8") as f:
                        f.write(migrated_text)
                    mig_reports.append((uploaded_files.name, mig_stats))

                # Sauvegarde de la position de cleavage
                if cleavage_start != 0:
                    config_path = os.path.join(entry_path, "_config.json")
                    with open(config_path, "w") as f:
                        json.dump(
                            {
                                "cleavage_start": int(cleavage_start),
                            },
                            f,
                            indent=2,
                        )

                # Save motifs as {entry_name}_motifs.csv
                if motif_file:
                    motifs_path = os.path.join(
                        entry_path, f"{entry_name.strip()}_motifs.csv"
                    )
                    with open(motifs_path, "wb") as f:
                        f.write(motif_file.getbuffer())

                st.session_state["form_key"] += 1
                st.session_state["success_msg"] = f"Entry '{entry_name}' added!"
                st.session_state["mig_reports"] = mig_reports
                st.rerun()

        if "success_msg" in st.session_state:
            st.success(st.session_state["success_msg"])
            del st.session_state["success_msg"]
        if "mig_reports" in st.session_state:
            for fname, mig_stats in st.session_state["mig_reports"]:
                with st.expander(f"Migration report - {fname}"):
                    st.write(
                        f"**{mig_stats['converted']} / {mig_stats['input']}** sequences kept"
                    )
                    st.write(
                        f"- {mig_stats['duplicates']} dropped - duplicate sequence"
                    )
                    st.write(f"- {mig_stats['no_genotype']} dropped - missing genotype")
                    st.write(
                        f"- {mig_stats['malformed']} dropped - malformed header (<7 fields)"
                    )
            del st.session_state["mig_reports"]
    st.divider()
    reset_button = st.button("Reset Cache", on_click=reset_app, type="primary")

    st.markdown(
        "##### **If you encounter any issues feel free to report them [here](https://github.com/Soupirr/NDV-genotyper/issues).**"
    )
# ============================================================================
# TITRE

if selection is None:
    st.title("Soupirr's Genotyper")
    st.markdown("""
    Analyze Virus sequences to identify genotypes and pathogenicity.
    """)
    st.divider()
# ============================================================================
# TABS

if not selection:
    st.info(
        "There is curently no reference sequences selected, use the sidebar panel to select one or add a new entry."
    )

else:
    st.title(f"{selection}")
    tabs = st.tabs(tab_labels)
    tab_analyze, tab_tree, tab_help, tab_val = tabs[0], tabs[1], tabs[2], tabs[3]
    selected_path = os.path.join(SEQ_FOLDER, selection)
    entry_config = load_entry_config(selected_path)
    with tab_help:
        stats_tab.render(selected_path, entry_config)
    with tab_analyze:
        analyze_tab.render(selected_path, entry_config)
    with tab_tree:
        tree_tab.render(selected_path)
    with tab_val:
        validation_tab.render(selected_path)
# ============================================================================
