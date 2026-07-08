"""Streamlit Web Application - GENOTYPER"""

import streamlit as st
from genotyper.tabs import tree_tab, analyze_tab, stats_tab, validation_tab
from genotyper.config import CUSTOM_CSS, SEQ_FOLDER
from genotyper.migration import migrate_fasta_text
import json
import os
import csv


# Page configuration
st.set_page_config(
    page_title="Soupirr's Genotyper",
    page_icon="misc/icon.png",
    layout="wide",
    initial_sidebar_state="expanded",
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

    motifs_files = [f for f in os.listdir(entry_path) if f.endswith("_motifs.csv")]
    if motifs_files:
        motifs_path = os.path.join(entry_path, motifs_files[0])
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


def reset_app():
    st.cache_resource.clear()
    st.cache_data.clear()
    st.session_state.pop("data_loaded", None)
    st.write("Cache cleared !")


if "sidebar_closed" not in st.session_state:
    st.session_state["sidebar_closed"] = False

if "sidebar_opened" not in st.session_state:
    st.session_state["sidebar_opened"] = True
    st.components.v1.html(
        """
        <script>
            setTimeout(function() {
                var sidebar = window.parent.document.querySelector('[data-testid="stSidebar"]');
                var collapseBtn = window.parent.document.querySelector('[data-testid="stSidebarCollapseButton"] button');
                if (sidebar && sidebar.getAttribute('aria-expanded') === 'false' && collapseBtn) {
                    collapseBtn.click();
                }
            }, 150);
        </script>
        """,
        height=0,
    )


with st.sidebar:
    st.link_button(
        "Documentation",
        "https://github.com/Soupirr/NDV-genotyper/blob/main/misc/QUICK_START.md",
    )
    entry = sorted(os.listdir(SEQ_FOLDER))
    st.divider()
    selection = st.selectbox(
        "Select an entry", entry, index=None, placeholder="Select an entry..."
    )

    if selection is not None and not st.session_state.get("sidebar_closed"):
        st.session_state["sidebar_closed"] = True
        st.components.v1.html(
            """
            <script>
                setTimeout(function() {
                    var btn = window.parent.document.querySelector('[data-testid="stSidebarCollapseButton"] button');
                    if (btn) btn.click();
                });
            </script>
            """,
            height=0,
        )

    if selection is None:
        st.session_state["sidebar_closed"] = False

    st.write("")
    with st.expander("Add new references datasets"):
        col1_side, col2_side = st.columns(2)

        with col1_side:
            if st.button("Mono"):
                st.session_state["mode"] = "mono"
                st.session_state.pop("gene_count", None)

        with col2_side:
            if st.button("Multi"):
                st.session_state["mode"] = "dual"
                st.session_state.pop("gene_count", None)

        st.divider()

        # ====================================================================================
        # =======================================MONO=========================================
        # ====================================================================================

        if st.session_state.get("mode") == "mono":
            entry_name = st.text_input(
                "Enter the name of the entry",
                key=f"name_{st.session_state['form_key']}",
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

        # ====================================================================================
        # =======================================DUAL=========================================
        # ====================================================================================

        if st.session_state.get("mode") == "dual":
            entry_name = st.text_input(
                "Enter the name of the entry",
                key=f"name_{st.session_state['form_key']}",
            )

            if "gene_count" not in st.session_state:
                st.session_state["gene_count"] = 1

            for i in range(st.session_state["gene_count"]):
                st.markdown(f"**Gene {i + 1}**")
                st.text_input("Gene name", key=f"gene_name_{i}")
                st.file_uploader(
                    "Fasta file(s)",
                    accept_multiple_files=True,
                    type=["fasta", "fas", "fa"],
                    key=f"gene_fasta_{i}",
                )
                st.text_input(
                    "Genotype pattern (ex: 'H\\d+' for H1/H2/H3/...)",
                    key=f"gene_patern_{i}",
                )
                st.write("")
                st.markdown(
                    "**Pathogenicity configuration** *(optional - leave blank to skip)*"
                )
                st.number_input(
                    "Cleavage start (optionnal)",
                    min_value=0,
                    step=1,
                    key=f"cleavage_{i}",
                )
                st.file_uploader(
                    "Motif file (optionnal)", type=["csv", "txt"], key=f"motif_{i}"
                )
                st.caption(
                    "Expected format: `motif,label,type` - e.g. `RRQKRF,VFcs-1,virulent`"
                )
                st.divider()

            if st.button("+ Add gene"):
                st.session_state["gene_count"] += 1
                st.rerun()

        # ====================================================================================
        # =====================================VALIDATION=====================================
        # ====================================================================================

        st.write("")

        adding_button = st.button("Add to the registry", type="primary")
        if adding_button:
            # =====================================================================MONO
            if st.session_state.get("mode") == "mono":
                if not upload_file:
                    st.error("Please select at least one FASTA file")
                elif not entry_name.strip():
                    st.error("Please enter a name")
                else:
                    entry_path = os.path.join(SEQ_FOLDER, entry_name.strip())
                    os.makedirs(entry_path, exist_ok=True)
                    mig_report = []
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
                        mig_report.append((uploaded_files.name, mig_stats))
                    if cleavage_start != 0:
                        config_path = os.path.join(entry_path, "_config.json")
                        with open(config_path, "w") as f:
                            json.dump(
                                {"cleavage_start": int(cleavage_start)}, f, indent=2
                            )
                    if motif_file:
                        motif_path = os.path.join(
                            entry_path, f"{entry_name.strip()}_motifs.csv"
                        )
                        with open(motif_path, "wb") as f:
                            f.write(motif_file.getbuffer())
                    st.session_state["form_key"] += 1
                    st.session_state["success_msg"] = f"Entry '{entry_name}' added !"
                    st.session_state["mig_reports"] = mig_report
                    st.rerun()

            # ======================================================================DUAL
            elif st.session_state.get("mode") == "dual":
                gene_count = st.session_state.get("gene_count", 1)
                entry_name_dual = entry_name.strip()

                # collecte des données de chaque gènes depuis la ss
                genes = []
                for i in range(gene_count):
                    genes.append(
                        {
                            "name": st.session_state.get(f"gene_name_{i}", "").strip(),
                            "files": st.session_state.get(f"gene_fasta_{i}", []),
                            "pattern": st.session_state.get(
                                f"gene_patern_{i}", ""
                            ).strip(),
                            "cleavage": st.session_state.get(f"cleavage_{i}", 0),
                            "motif": st.session_state.get(f"motif_{i}", None),
                        }
                    )

                # Validation des données
                if not entry_name_dual:
                    st.error("Please enter an entry name")
                elif any(not g["name"] for g in genes):
                    st.error("Please enter a name for each genes")
                elif any(not g["files"] for g in genes):
                    st.error("Please upload at least one FASTA per gene")

                else:
                    entry_path = os.path.join(SEQ_FOLDER, entry_name_dual)
                    os.makedirs(entry_path, exist_ok=True)
                    mig_report = []

                    for g in genes:
                        gene_path = os.path.join(entry_path, g["name"])
                        os.makedirs(gene_path, exist_ok=True)

                        # FASTA
                        for uploaded_files in g["files"]:
                            dest_path = os.path.join(gene_path, uploaded_files.name)
                            raw_text = (
                                uploaded_files.getbuffer()
                                .tobytes()
                                .decode("utf-8-sig", errors="ignore")
                            )
                            migrated_text, mig_stats = migrate_fasta_text(raw_text)
                            with open(dest_path, "w", encoding="utf-8") as f:
                                f.write(migrated_text)
                            mig_report.append(
                                (f"{g['name']}/{uploaded_files.name}", mig_stats)
                            )

                        # Config json
                        config = {}
                        if g["cleavage"] != 0:
                            config["cleavage_start"] = int(g["cleavage"])
                        if g["pattern"]:
                            config["genotype_pattern"] = g["pattern"]
                        if config:
                            with open(
                                os.path.join(gene_path, "_config.json"), "w"
                            ) as f:
                                json.dump(config, f, indent=2)

                        # Fichier motif
                        if g["motif"]:
                            motif_path = os.path.join(
                                gene_path, f"{g['name']}_motifs.csv"
                            )
                            with open(motif_path, "wb") as f:
                                f.write(g["motif"].getbuffer())

                    st.session_state["form_key"] += 1
                    st.session_state["success_msg"] = (
                        f"Entry '{entry_name_dual}' added!"
                    )
                    st.session_state["mig_reports"] = mig_report
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
