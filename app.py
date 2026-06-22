"""
Newcastle Disease Virus F Gene Genotyper
Streamlit Web Application
"""

import streamlit as st
import pandas as pd
import time
import plotly.graph_objects as go
from ndv_genotyper.analyzer import (
    FASTAParser,
    analyze_newcastle_sequence,
    unpack_top_match,
    get_class,
    SequenceSimilarity,
)
from ndv_genotyper.tabs import map_tab, tree_tab
import os
import platform

if platform.system() == "Windows":
    import winsound

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "data")
SEQ_FOLDER = os.path.join(DATA_FOLDER, "sequences")
LOCATION_FOLDER = os.path.join(DATA_FOLDER, "locations")
PATHO_FOLDER = os.path.join(DATA_FOLDER, "pathogenicity")

# palette importé depuis ColorBrewer (Spectral 11 + Set1/Dark2 extensions)
# les couleurs pâles du centre Spectral (#fee08b, #ffffbf, #e6f598, #abdda4)
# ont été remplacées par des équivalents saturés — Spectral est conçu pour fond blanc
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


st.markdown(
    """
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
""",
    unsafe_allow_html=True,
)


# ============================================================================
# MAIN TABS
# ============================================================================

# def des onglets
tab_analyze, tab_tree, tab_map, help_tab = st.tabs(
    [
        "Analyze Sequences",
        "Phylogenetic Trees",
        "Global Distribution Map",
        "Help - Statistics",
    ]
)


# ============================================================================
# TAB 1: ANALYZE SEQUENCE
# ============================================================================

with tab_analyze:
    # Charge les datasets de référence / reset les dictionnaires
    @st.cache_resource
    def load_all_references():
        combined_sequences = {}
        errors = []

        # Charge automatiquement tous les fichiers .fas du dossier "sequences"
        fas_files = [f for f in os.listdir(SEQ_FOLDER) if f.endswith(".fas")]

        for filename in fas_files:
            try:
                seqs = FASTAParser.parse_file(os.path.join(SEQ_FOLDER, filename))
                if seqs:
                    combined_sequences.update(seqs)
            except Exception as e:
                errors.append(f"{filename}: {e}")

        return combined_sequences, len(fas_files), len(combined_sequences), errors

    references, files_count, total_count, load_errors = load_all_references()

    st.header(
        f"Sequence Analyzer Tool \u00a0\u00a0\u00a0({files_count} Databases Loaded)",
        help="files can be found in '.../data/sequences/*'",
    )

    if (
        not references or len(references) == 0
    ):  # Failsafe si il n'y a pas les fasta dans le dossier
        st.error("Error: Could not load any reference datasets")
        for error in load_errors:
            st.warning(f"⚠ {error}")
        st.info("""
        **How to fix:**
        1. Download the NDV reference datasets from: https://github.com/NDVconsortium/NDV_Sequence_Datasets
        2. Place these files in the "sequences" folder
        3. Restart the app
        """)
        st.stop()

    st.divider()

    # layout: Configuration sur la gauche, Input/Results sur la droite
    col_config, col_content = st.columns([0.2, 0.8])

    with col_config:
        st.subheader("⚙ Configuration")

        similarity_method = st.radio(  # bouton radio pour selectionner la méthode d'analyse
            "Similarity Method",
            ["Hamming (fast)", "Pairwise (accurate)"],
            help="Hamming: simple mismatch count. Pairwise: accounts for insertions/deletions",
        )

        # Alerte pour la méthode de Leveshtein
        if "Pairwise" in similarity_method:
            st.warning(
                "⚠ **Pairwise Method**\n\n"
                "This method is more accurate but slower. "
                "Analysis can take **5-30+ minutes** depending on number of sequences and computer power. "
                "Use **Hamming** for faster results."
            )

        st.write("")

        # Scroll pour choisir le nombre de Matches à afficher
        top_n_matches = st.slider(
            "Number of Top Matches",
            min_value=1,
            max_value=10,
            value=5,
            help="Show this many best-matching genotypes  (Default = 5)",
        )

        st.write("")

        show_matrix = st.checkbox(
            "Compute Similarity Matrix",
            value=True,
            help="This should be turned off when analysing more than ~20 sequences at once",
        )

    with col_content:
        st.subheader("Input Sequence")

        input_fasta = None

        def clear_input():
            st.session_state["fasta_input_area"] = ""
            st.session_state.pop("all_sequences", None)  # clear sequences
            st.session_state.pop("all_results", None)  # clear results
            st.session_state.pop("trees", None)  # clear cached trees
            st.session_state.pop("load_tree", None)
            st.session_state.pop("n_neighbours", None)

        with st.form("analyze_form", border=False):
            input_tab1, input_tab2 = st.tabs(["Paste FASTA", "Upload File"])

            # Input sous forme de texte
            with input_tab1:
                input_fasta = st.text_area(
                    "Paste your Newcastle Disease F gene sequence in FASTA format:",
                    height=200,
                    placeholder=">your_sequence_name\nATGGGCTCCAGATCCTCTAC...",
                    help="Should be a complete F gene sequence (~1662 bp)",
                    key="fasta_input_area",
                )

            # Input via upload de fichier
            with input_tab2:
                uploaded_file = st.file_uploader(  # fichier ouvert en mode binaire.
                    "Or upload a FASTA file", type=["fasta", "fas", "fa", "txt"]
                )
                if uploaded_file:
                    input_fasta = uploaded_file.read().decode(
                        "utf-8"
                    )  # transforme le fichier uploadé (bytes) en texte

            # bouton d'analyse
            col_btn1, col_btn2 = st.columns([0.7, 0.3])
            with col_btn1:
                analyze_button = st.form_submit_button(
                    "Analyze Sequences", type="primary", width="stretch"
                )
            # bouton de reset
            with col_btn2:
                st.form_submit_button(
                    "Clear Input",
                    width="stretch",
                    on_click=clear_input,
                )

        if analyze_button:
            st.session_state.pop("load_tree", None)
            st.session_state.pop("n_neighbours", None)
            if not input_fasta or input_fasta.strip() == "":
                st.error("Please provide a sequence to analyze")
            elif not input_fasta.strip().startswith(">"):
                st.error("Please provide a correct sequence")
            else:
                all_sequences = FASTAParser.parse_text(input_fasta)
                all_results = []
                start_time = time.time()
                method = "hamming" if "Hamming" in similarity_method else "pairwise"

                # ajout d'une barre de progression
                progress_bar = st.progress(0, text="Starting analysis...")
                all_results = []
                for index, (header, sequence) in enumerate(all_sequences.items()):
                    progress = (index + 1) / len(all_sequences)
                    progress_bar.progress(
                        progress,
                        text=f"Analyzing {index + 1}/{len(all_sequences)}: {header[:30]}...",
                    )

                    single_fasta = f">{header}\n{sequence}"
                    result = analyze_newcastle_sequence(
                        input_fasta=single_fasta,
                        reference_sequences=references,
                        top_matches=top_n_matches,
                        similarity_method=method,
                    )
                    all_results.append(result)
                if platform.system() == "Windows":
                    winsound.PlaySound(
                        "misc/notification.wav",
                        winsound.SND_FILENAME | winsound.SND_ASYNC,
                    )
                progress_bar.progress(1, text="Analysis complete!")
                progress_bar.empty()

                # sauvegarde des informations de l'analyse
                st.session_state["all_results"] = all_results
                st.session_state["all_sequences"] = all_sequences
                st.session_state["elapsed_time"] = time.time() - start_time
                st.session_state["method"] = method

        if "all_results" in st.session_state and st.session_state["all_results"]:
            all_results = st.session_state["all_results"]
            all_sequences = st.session_state.get("all_sequences", {})
            method = st.session_state.get("method", "hamming")
            elapsed_time = st.session_state.get("elapsed_time", 0)

            st.success(f"Analysis of {len(all_results)} sequence(s) completed!")
            st.info(f"Analysis took **{elapsed_time:.2f} seconds**")

            # Display des résultats
            seq_tabs = st.tabs([r["input_header"][:30] for r in all_results])

            for i, (tab, results) in enumerate(zip(seq_tabs, all_results)):
                with tab:
                    st.subheader("Sequence Information")
                    col1, col2 = st.columns([0.8, 0.2])  # 1er colonne (80%)
                    with col1:
                        st.write(f"**Header:** {results['input_header']}")  # header
                    with col2:
                        st.write(f"**Length:** {results['sequence_length']} bp")  # bp

                    st.divider()

                    st.subheader("Genotype Identification")
                    matches = results["genotype_matches"]
                    # import les matches de results

                    if matches:
                        top_match_raw = matches[0]
                        top_match = unpack_top_match(top_match_raw)
                        # on trie les matches par ordre croissant de similutude donc le top matches est [0]
                        col0, col1, col2, col3, col4 = st.columns(5)
                        with col0:
                            st.metric("Class", get_class(top_match["genotype"]))
                        with col1:
                            st.metric(
                                "Best Match Genotype",
                                top_match["genotype"],
                                f"{top_match['avg_similarity']}%",
                            )
                        with col2:
                            st.metric(
                                "Reference Sequences",
                                top_match["sample_count"],
                            )
                        with col3:
                            # affiche le header du top match
                            st.metric(
                                "Average Match Score",
                                f"{top_match['avg_similarity']}%",
                                "_".join(
                                    top_match["best_header"]
                                    .replace(">", "")
                                    .split("_")[2:6]
                                ),
                            )
                        with col4:
                            st.metric(
                                "Similarity Method",
                                "Hamming"
                                if "Hamming" in similarity_method
                                else "Pairwise",
                            )

                        if (
                            method == "hamming"
                            and matches
                            and top_match["avg_similarity"] < 95
                        ):
                            st.warning(
                                "⚠ Average Similarity score is below 95% — result may not be interpretable with the Hamming method. "
                                "Consider using the Pairwise method for more accurate results."
                            )

                        # Tableau avec les Top Matches
                        st.write("**Top Matching Genotypes:**")
                        # Transformation en DataFrame
                        matches_df = pd.DataFrame(
                            matches,
                            columns=[
                                "Genotype",
                                "Avg Similarity (%)",
                                "Reference Count",
                                "Best Match Header",
                                "Best Match Score (%)",
                            ],
                        )
                        matches_df = matches_df[  # Reorganisation des colonnes
                            [
                                "Genotype",
                                "Avg Similarity (%)",
                                "Best Match Score (%)",
                                "Reference Count",
                            ]
                        ]
                        st.dataframe(
                            matches_df,
                            width="stretch",
                            hide_index=True,
                        )
                        st.caption(
                            "Reference Count = number of reference sequences for each genotype in reference Dataset"
                        )

                    st.divider()

                    st.subheader("Pathogenicity Analysis")

                    cadre_1, cadre_2, cadre_3 = st.tabs(
                        [
                            "Main Reading Frame",
                            "Reading Frame +1",
                            "Reading Frame -1",
                        ]
                    )

                    def display_reading_frame(cleavage: dict, frame_label: str):
                        if cleavage["cleavage_region_found"]:
                            motif = (
                                cleavage["motif_type"]
                                if cleavage["motif_type"]
                                and len(cleavage["motif_type"]) == 6
                                else cleavage["cleavage_protein"][:8]
                            )
                            protein = cleavage["cleavage_protein"]

                            col1, col2, col3 = st.columns([2, 1, 1])

                            with col1:
                                st.write("**Cleavage Region (protein):**")
                                if motif and len(motif) == 6 and motif in protein:
                                    idx = protein.index(motif)
                                    before = protein[:idx]
                                    after = protein[idx + len(motif) :]
                                    color = (
                                        "#FF4444"
                                        if "Virulent" in cleavage["pathogenicity"]
                                        else "#2ECC71"
                                    )
                                    st.markdown(
                                        f"<div style='font-size:20px; text-align:center; font-family:monospace; background-color:#0e1117; padding:8px; border-radius:4px; border:1px solid rgba(255,255,255,0.2);'>{before}<span style='color:{color}'>{motif}</span>{after}</div>",
                                        unsafe_allow_html=True,
                                    )
                                else:
                                    # si pas de motifs trouvés
                                    st.markdown(
                                        f"<div style='font-size:20px; text-align:center; font-family:monospace; background-color:#0e1117; padding:8px; border-radius:4px; border:1px solid rgba(255,255,255,0.2); color:#868d2f;'>{protein}</div>",
                                        unsafe_allow_html=True,
                                    )
                            with col2:
                                st.write("**Motif Found:**")
                                if (
                                    cleavage["motif_type"]
                                    and len(cleavage["motif_type"]) == 6
                                ):
                                    st.code(cleavage["motif_type"], language="text")
                                else:
                                    st.code("Not found", language="text")
                            with col3:
                                st.write("**Motif Category:**")
                                st.code(
                                    f"{cleavage['motif_category'] or 'Unknown'}",
                                    language="text",
                                )

                            pathogenicity = cleavage["pathogenicity"]
                            confidence_path = cleavage["confidence"]

                            # Ici on utilise error (et success) pour afficher la zone de texte en rouge ou vert (purement cosmétique)
                            if "Virulent" in pathogenicity:
                                st.error(f"Virulent (Confidence: {confidence_path})")
                                st.write(
                                    "This sequence shows characteristics of a virulent strain with polybasic cleavage site (≥3 basic residues at positions 113-116 + F at position 117)."
                                )
                            elif "Low-virulence" in pathogenicity:
                                st.success(
                                    f"Low-virulence (Confidence: {confidence_path})"
                                )
                                st.write(
                                    "This sequence shows characteristics of a low-virulence strain with monobasic or L117 cleavage site."
                                )
                            else:
                                st.warning(
                                    f"Undetermined (Confidence: {confidence_path})"
                                )
                                st.write(
                                    "No known cleavage site motif was found in the expected region. "
                                    "The sequence may be incomplete, contain too many indels, or represent an unusual strain."
                                )
                        else:
                            st.warning(
                                "Could not analyze cleavage site. Sequence may be incomplete."
                            )

                    with cadre_1:
                        display_reading_frame(results["cleavage_main"], "Main")
                    with cadre_2:
                        display_reading_frame(results["cleavage_plus_one"], "+1")
                    with cadre_3:
                        display_reading_frame(results["cleavage_minus_one"], "-1")

                    st.caption(
                        "Pathogenicity prediction is based on cleavage site motif matching (Wang et al., 2017) and should be interpreted as indicative only. Results must be validated by biological assays before any clinical or regulatory conclusion."
                    )

            # Matrice de similitude
            if len(all_sequences) > 1 and show_matrix:
                st.divider()
                st.subheader("Input Sequences Comparison Matrix")

                seq_headers_short = list(all_sequences.keys())
                seq_headers_full = list(all_sequences.keys())
                matrix = []

                for header1, seq1 in all_sequences.items():
                    row = []
                    for header2, seq2 in all_sequences.items():
                        similarity = SequenceSimilarity.pairwise_similarity(seq1, seq2)
                        row.append(similarity)
                    matrix.append(row)

                hover_text = [
                    [
                        f"{seq_headers_full[i]}<br>{seq_headers_full[j]}<br>Similarity: {matrix[i][j]:.1f}%"
                        for j in range(len(seq_headers_full))
                    ]
                    for i in range(len(seq_headers_full))
                ]

                fig_matrix = go.Figure(
                    data=go.Heatmap(
                        z=matrix,
                        x=[h[:10] for h in seq_headers_short],
                        y=[h[:10] for h in seq_headers_short],
                        colorscale="RdYlGn_r",
                        zmid=85,
                        text=[[f"{val:.2f}%" for val in row] for row in matrix],
                        texttemplate="%{text}",
                        textfont={"size": 10},
                        colorbar=dict(title="Similarity %"),
                        hovertext=hover_text,
                        hovertemplate="%{hovertext}<extra></extra>",
                    )
                )
                fig_matrix.update_layout(
                    height=730,
                    width=450,
                    plot_bgcolor="rgba(0,0,0,0)",
                    paper_bgcolor="rgba(0,0,0,0)",
                    margin=dict(t=40, l=50, r=50, b=50),
                )
                st.plotly_chart(fig_matrix, width="stretch")

            st.divider()
            st.subheader("Full Analysis Details:")

            export_rows = []  # Tableau de fin
            for results in all_results:
                matches = results["genotype_matches"]
                frames = [
                    (results["cleavage_main"], "Main"),
                    (results["cleavage_plus_one"], "+1"),
                    (results["cleavage_minus_one"], "-1"),
                ]
                best_cleavage = None
                best_frame_label = "N/A"

                for cleavage, label in frames:
                    if cleavage["pathogenicity"] != "Undetermined":
                        best_cleavage = cleavage
                        best_frame_label = label
                        break

                export_rows.append(
                    {
                        "Header": results["input_header"],
                        "Length": results["sequence_length"],
                        "Class": get_class(matches[0][0]) if matches else "N/A",
                        "Best Genotype": matches[0][0] if matches else "N/A",
                        "Similarity %": matches[0][1] if matches else "N/A",
                        "Pathogenicity": best_cleavage["pathogenicity"]
                        if best_cleavage
                        else "Undetermined",
                        "Cleavage Motif": f"{best_cleavage['motif_type']}"
                        if best_cleavage
                        else "N/A",
                        "Motif Class": f"{best_cleavage['motif_category']}"
                        if best_cleavage
                        else "N/A",
                        "Reading Frame": best_frame_label
                        if best_frame_label
                        else "N/A",
                    }
                )

            export_df = pd.DataFrame(export_rows)
            st.dataframe(export_df, width="stretch", hide_index=True)


# ============================================================================
# TAB 2: PHYLOGENETIC TREES
# ============================================================================

with tab_tree:
    tree_tab.render(references)

# ============================================================================
# TAB 3: GLOBAL DISTRIBUTION MAP
# ============================================================================

with tab_map:
    map_tab.render()

# ============================================================================
# TAB 4: Stats/Help
# ============================================================================

with help_tab:
    st.markdown(
        "##### **If you encounter any issues feel free to report them [here](https://github.com/Soupirr/NDV-genotyper/issues).**"
    )
    st.divider()

    (
        stat_tab,
        info_tab,
        path_tab,
    ) = st.tabs(["Statistics", "Informations", "Pathogenicity"])

    with info_tab:
        with open("misc/QUICK_START.md", "r", encoding="utf-8") as f:
            informations = f.read()

        st.markdown(informations)

    with path_tab:
        with open("misc/PATHOGENICITY_CRITERIA.md", "r", encoding="utf-8") as p:
            pathogénicité = p.read()

        st.markdown(pathogénicité)

    with stat_tab:
        st.header("Database Statistics")

        db_references, db_files_count, db_total_count, _ = load_all_references()

        # analyse de tout les champs de chaque header
        genotype_counts = {}
        class_counts = {"Class I": 0, "Class II": 0}
        host_counts = {}
        year_counts = {}

        # Normalisation des noms d'hôtes
        # Pour ajouter un hôte rajouter une ligne raw_host,normalized_host  dans le CSV
        _df_host_map = pd.read_csv(
            os.path.join(DATA_FOLDER, "hosts", "host_normalize.csv")
        )
        HOST_NORMALIZE = dict(
            zip(_df_host_map["raw_host"].str.lower(), _df_host_map["normalized_host"])
        )
        df_patho = pd.read_csv(
            os.path.join(PATHO_FOLDER, "NDV_F_class_I_and_II_2022_Pathogenicity.csv")
        )

        KNOWN_COUNTRIES = set(map_tab.df_world["country"].str.lower()) | {
            "china",
            "usa",
            "russia",
        }

        # US abbreviations
        US_STATE_ABBREVS = set(map_tab.df_us["abbreviation_us"].str.lower())

        for header in db_references.keys():
            parts = header.split("_")
            if len(parts) < 2:
                continue

            # Filtrage Genotypes et Class
            genotype = parts[1]
            genotype_counts[genotype] = genotype_counts.get(genotype, 0) + 1
            if genotype.startswith("1"):
                class_counts["Class I"] += 1
            else:
                class_counts["Class II"] += 1

            # Filtrage des animaux
            host_idx = 5
            if len(parts) > 5 and parts[5] == "1":
                host_idx = 6  # parfoit un "1" inexpliqué dans le header
            if len(parts) > host_idx:
                raw_host = parts[host_idx].lower()
                if (
                    not raw_host
                    or raw_host in KNOWN_COUNTRIES
                    or raw_host in US_STATE_ABBREVS
                    or raw_host.startswith("apmv")
                    or raw_host.startswith("ndv")
                    or raw_host.startswith("aoav")
                ):
                    # si hôte non reconnus
                    host = "Unspecified"
                else:
                    host = HOST_NORMALIZE.get(raw_host, "Other")
                host_counts[host] = host_counts.get(host, 0) + 1

            # Filtrage des années
            last = parts[-1]
            if last.isdigit() and len(last) == 4:
                year_counts[last] = year_counts.get(last, 0) + 1

        df_genotypes = pd.DataFrame(
            sorted(genotype_counts.items(), key=lambda x: x[1], reverse=True),
            columns=["Genotype", "Sequences"],
        )

        _df_map_stats = map_tab.build_map_dataframe()
        _df_map_stats = _df_map_stats[_df_map_stats["label"] != "Unknown"].copy()
        _df_map_stats["country"] = _df_map_stats["label"].apply(
            lambda x: x.split(",")[-1].strip()
        )
        country_counts = _df_map_stats["country"].value_counts().to_dict()

        # Statistique de bases
        years_with_data = [int(y) for y in year_counts if year_counts[y] > 0]
        col2, col3, col4, col5, col6 = st.columns(5)
        col2.metric("Total Sequences", db_total_count)
        col3.metric("Unique Sub-Genotypes", len(genotype_counts))
        col4.metric("Class I Sequences", class_counts["Class I"])
        col5.metric("Class II Sequences", class_counts["Class II"])
        col6.metric(
            "Year Range",
            f"{min(years_with_data)}–{max(years_with_data)}"
            if years_with_data
            else "N/A",
        )

        st.divider()

        # sequences par genotype et class pie chart
        col_bar, col_pie = st.columns([3, 1])

        with col_bar:
            st.subheader("Sequences per Genotype")
            fig_bar = go.Figure(
                go.Bar(
                    x=df_genotypes["Genotype"],
                    y=df_genotypes["Sequences"],
                    marker=dict(
                        color=df_genotypes["Sequences"],
                        colorscale=[[0, "#0099cc"], [1, "#00c9a7"]],
                        showscale=False,
                    ),
                    text=df_genotypes["Sequences"],
                    textposition="outside",
                )
            )
            fig_bar.update_layout(
                xaxis_title="Genotype",
                yaxis_title="Nb. of Sequences",
                height=420,
                margin=dict(l=0, r=0, t=0, b=0),
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                font=dict(color="white"),
                xaxis=dict(tickangle=-45),
            )
            st.plotly_chart(fig_bar, width="stretch")

        with col_pie:
            st.subheader("Class Distribution")
            fig_pie = go.Figure(
                go.Pie(
                    labels=list(class_counts.keys()),
                    values=list(class_counts.values()),
                    hole=0.45,
                    marker=dict(colors=["#0099cc", "#00c9a7"]),
                    textinfo="label+percent",
                    textfont=dict(color="white"),
                )
            )
            fig_pie.update_layout(
                height=420,
                margin=dict(l=0, r=0, t=0, b=0),
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                font=dict(color="white"),
                showlegend=False,
            )
            st.plotly_chart(fig_pie, width="stretch")

        st.divider()

        # distribution temporel des séquences
        st.subheader("Sequences by Year")
        df_years = pd.DataFrame(
            sorted(
                [
                    (int(y), c)
                    for y, c in year_counts.items()
                    if y.isdigit() and len(y) == 4
                ]
            ),  # code assez restrictif qui ne prend que un nombre si le header finis par 4 nombres (c'est la galère)
            columns=["Year", "Sequences"],
        )
        fig_year = go.Figure(
            go.Bar(
                x=df_years["Year"],
                y=df_years["Sequences"],
                marker=dict(color="#00c9a7"),
                text=df_years["Sequences"],
                textposition="outside",
            )
        )
        fig_year.update_layout(
            xaxis_title="Year",
            yaxis_title="Nb. of Sequences",
            height=350,
            margin=dict(l=0, r=0, t=0, b=0),
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            font=dict(color="white"),
            xaxis=dict(dtick=2),
        )
        st.plotly_chart(fig_year, width="stretch")

        st.divider()

        # pie chart des hôtes et graph des pays
        col_host, col_country = st.columns([4, 6])

        with col_host:
            st.subheader("Host Distribution")
            host_items = sorted(host_counts.items(), key=lambda x: x[1], reverse=True)
            host_items = [(k, v) for k, v in host_items if k == "Unspecified"] + [
                (k, v) for k, v in host_items if k != "Unspecified"
            ]
            df_hosts = pd.DataFrame(host_items, columns=["Host", "Sequences"])

            # PALETTE définie en haut du fichier (ColorBrewer)
            pie_colors = []
            palette_i = 0
            for host in df_hosts["Host"]:
                if host == "Unspecified":
                    pie_colors.append("#888888")
                elif host == "Other":
                    pie_colors.append("#555555")
                else:
                    pie_colors.append(PALETTE[palette_i % len(PALETTE)])
                    palette_i += 1

            fig_host = go.Figure(
                go.Pie(
                    labels=df_hosts["Host"],
                    values=df_hosts["Sequences"],
                    hole=0.4,
                    marker=dict(colors=pie_colors),
                    textinfo="label+percent",
                    textposition="inside",
                    textfont=dict(color="white"),
                )
            )
            fig_host.update_layout(
                height=380,
                margin=dict(l=0, r=0, t=30, b=0),
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                font=dict(color="white"),
                showlegend=False,
            )
            st.plotly_chart(fig_host, width="stretch")

        # affichage du graph des pays
        with col_country:
            st.subheader("Top 15 Countries")
            df_countries = pd.DataFrame(
                sorted(country_counts.items(), key=lambda x: x[1], reverse=True)[:15],
                columns=["Country", "Sequences"],
            )
            fig_country = go.Figure(
                go.Bar(
                    x=df_countries["Sequences"],
                    y=df_countries["Country"],
                    orientation="h",
                    marker=dict(
                        color=df_countries["Sequences"],
                        colorscale=[[0, "#00b0ea"], [1, "#00b194"]],
                        showscale=False,
                    ),
                    text=df_countries["Sequences"],
                    textposition="outside",
                )
            )
            fig_country.update_layout(
                xaxis_title="Number of Sequences",
                height=380,
                margin=dict(l=0, r=30, t=30, b=0),
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                font=dict(color="white"),
                yaxis=dict(autorange="reversed"),
            )
            st.plotly_chart(fig_country, width="stretch")

        st.divider()

        # check-up de la santé actuel du jeu de donnée
        # trop peu de séquence dans un génotype n'est pas un bonne chose pour les analyses
        # pour ça on créer un dataframe avec les génotypes et leurs comptes
        st.subheader("⚠ Database Coverage Health")
        THRESHOLD = 15
        df_health = df_genotypes.copy()
        df_health["Status"] = df_health["Sequences"].apply(
            lambda n: "Good" if n >= THRESHOLD else ("Low" if n >= 5 else "Critical")
        )
        df_health["Class"] = df_health["Genotype"].apply(
            lambda g: "Class I" if g.startswith("1") else "Class II"
        )
        df_health["% of Total"] = (df_health["Sequences"] / db_total_count * 100).map(
            "{:.1f}%".format
        )
        df_health = df_health[
            ["Genotype", "Class", "Sequences", "% of Total", "Status"]
        ]

        critical = df_health[df_health["Status"] == "Critical"]
        low = df_health[df_health["Status"] == "Low"]

        c1, c2, c3 = st.columns(3)
        c1.metric(
            "Well covered (≥15 seq)",
            len(df_health[df_health["Status"] == "Good"]),
        )
        c2.metric("Low coverage (5–14 seq)", len(low))
        c3.metric("Critical (<5 seq)", len(critical))

        if not critical.empty:
            st.markdown("**Genotypes needing urgent attention:**")
            st.dataframe(critical, width="stretch", hide_index=True)

        # affichage du dataframe si nécessaire
        with st.expander("View full coverage table"):
            st.dataframe(df_health, width="stretch", hide_index=True)

        st.divider()

        st.subheader("Genotypes Pathogenicity")

        fig_gen_patho = go.Figure()
        patho_counts = (
            df_patho.groupby(["Best Genotype", "Pathogenicity"])
            .size()
            .reset_index(name="count")
        )
        for patho_type in patho_counts["Pathogenicity"].unique():
            df_p = patho_counts[patho_counts["Pathogenicity"] == patho_type]
            fig_gen_patho.add_trace(
                go.Bar(
                    name=patho_type,
                    x=df_p["Best Genotype"],
                    y=df_p["count"],
                    hovertemplate="%{x}<br>%{y:.2f}%<br>%{fullData.name}<extra></extra>",
                    marker_color={
                        "Virulent": "#FC3C3C",
                        "Low-virulence": "#2ECC71",
                        "Undetermined": "#888888",
                    }.get(patho_type, "#888888"),
                )
            )
        fig_gen_patho.update_layout(
            barmode="stack",
            barnorm="percent",
            legend=dict(yanchor="bottom", y=-0.5, xanchor="left", x=0.01),
            xaxis_title="Genotypes",
            yaxis_title="% Virulence",
            height=350,
            margin=dict(l=0, r=0, t=30, b=0),
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            hoverlabel=dict(
                bgcolor="#0c1a24",
                font_color="white",
                bordercolor="rgba(255,255,255,0.2)",
            ),
        )
        st.plotly_chart(fig_gen_patho, width="stretch")

        st.divider()

        col_patho_1, col_patho_2 = st.columns([2, 1])

        with col_patho_1:
            st.subheader("Motifs Distribution")
            motif_counts = (
                df_patho.groupby(["Cleavage Motif", "Motif Class"])
                .size()
                .reset_index(name="Count")
            )
            motif_counts = motif_counts[
                motif_counts["Cleavage Motif"] != "N/A"
            ].sort_values("Count", ascending=False)
            st.dataframe(motif_counts, hide_index=True, width="stretch")

        with col_patho_2:
            st.subheader("Virulence Distribution")
            vir_counts = df_patho["Pathogenicity"].value_counts()
            fig_vir_pie = go.Figure(
                go.Pie(
                    labels=vir_counts.index,
                    values=vir_counts.values,
                    hole=0.45,
                    textposition="inside",
                    marker=dict(
                        colors=[
                            {
                                "Virulent": "#FF4444",
                                "Low-virulence": "#2ECC71",
                                "Undetermined": "#888888",
                            }.get(lll, "#888888")
                            for lll in vir_counts.index
                        ]
                    ),
                    textinfo="label+percent",
                    textfont=dict(color="white"),
                )
            )
            fig_vir_pie.update_layout(
                height=350,
                margin=dict(l=0, r=0, t=0, b=0),
                paper_bgcolor="rgba(0,0,0,0)",
                font=dict(color="white"),
                showlegend=False,
            )
            st.plotly_chart(fig_vir_pie, width="stretch")
