"""
Newcastle Disease Virus F Gene Genotyper
Streamlit Web Application
"""

import streamlit as st
import winsound
import pandas as pd
import time
import plotly.graph_objects as go
from analyzer import (
    FASTAParser,
    analyze_newcastle_sequence,
    unpack_top_match,
    get_class,
)
from analyzer import (
    build_tree_fasttree,
    get_color,
    align_sequences_mafft,
    find_closest_neighbours,
    write_temp_fasta,
)
import os
import pydeck as pdk

DATA_FOLDER = os.path.join(os.path.dirname(__file__), "data")
SEQ_FOLDER = os.path.join(DATA_FOLDER, "sequences")
LOCATION_FOLDER = os.path.join(DATA_FOLDER, "locations")

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
    st.header("Sequence Analyzer Tool")

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

    # Count unique genotypes
    unique_genotypes = (
        set()
    )  # set() permet d'ajouter des données en ignorant les doublons
    for header in references.keys():  # parcours les headers du dict du jeu de données
        parts = header.split("_")  # coupe au niveau des "_"
        if len(parts) >= 2:  # vérifie qu'il y a au moins 2 morceaux après le découpage
            genotype = parts[1]  # prend le 2ème morceau (le génotype)
            unique_genotypes.add(genotype)  # l'ajoute à "unique_genotype"

    # Display dataset info
    col1, col2, _ = st.columns([2, 2, 6])
    with col1:
        st.metric(
            "Databases Loaded",
            files_count,
            help="files can be found in '.../data/sequences/*'",
        )
    with col2:
        st.metric(
            "Total Sequences",
            len(references),
            help="unique genotypes names can be found in '.../data/genotypes.txt'",
        )

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

        # Scroll pour choisir le nombre de Matches à afficher
        top_n_matches = st.slider(
            "Number of Top Matches",
            min_value=1,
            max_value=10,
            value=3,
            help="Show this many best-matching genotypes  (Default = 3)",
        )

    with col_content:
        st.subheader("Input Sequence")

        input_tab1, input_tab2 = st.tabs(["Paste FASTA", "Upload File"])

        input_fasta = None

        def clear_input():
            st.session_state["fasta_input_area"] = ""
            st.session_state.pop("all_sequences", None)  # clear sequences
            st.session_state.pop("all_results", None)  # clear results
            st.session_state.pop("trees", None)  # clear cached trees

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
            analyze_button = st.button(
                "Analyze Sequences", type="primary", width="stretch"
            )
        # bouton de reset
        with col_btn2:
            st.button(
                "Clear Input",
                width="stretch",
                on_click=clear_input,
            )

        if analyze_button:
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
                winsound.PlaySound(
                    "misc/notification.wav", winsound.SND_FILENAME | winsound.SND_ASYNC
                )
                progress_bar.progress(1, text="Analysis complete!")

                # sauvegarde des informations de l'analyse
                st.session_state["all_results"] = all_results
                st.session_state["all_sequences"] = all_sequences
                st.session_state["elapsed_time"] = time.time() - start_time
                st.session_state["method"] = method

        if "all_results" in st.session_state and st.session_state["all_results"]:
            all_results = st.session_state["all_results"]
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

                        st.divider()

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

                    # import les résultats de l'analyse
                    cleavage = results["cleavage_analysis"]

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
                        confidence = cleavage["confidence"]

                        # Ici on utilise error (et success) pour afficher la zone de texte en rouge ou vert (purement cosmétique)
                        if "Likely Virulent" in pathogenicity:
                            st.error(f"Likely Virulent (Confidence: {confidence})")
                            st.write(
                                "This sequence shows characteristics of a virulent strain with polybasic cleavage site (≥3 basic residues at positions 113-116 + F at position 117)."
                            )
                        elif "Likely Low-virulence" in pathogenicity:
                            st.success(
                                f"Likely Low-virulence (Confidence: {confidence})"
                            )
                            st.write(
                                "This sequence shows characteristics of a low-virulence strain with monobasic or L117 cleavage site."
                            )
                        else:
                            st.warning(f"Undetermined (Confidence: {confidence})")
                            st.write(
                                "No known cleavage site motif was found in the expected region. "
                                "The sequence may be incomplete, contain too many indels, or represent an unusual strain."
                            )
                    else:
                        st.warning(
                            "Could not analyze cleavage site. Sequence may be incomplete."
                        )

            # Exportation des résultats dans un fichier .csv
            st.divider()
            st.subheader("Full Analysis Details:")

            export_rows = []  # Tableau de fin
            for results in all_results:
                matches = results["genotype_matches"]
                cleavage = results["cleavage_analysis"]
                export_rows.append(
                    {
                        "Input Sequence": results["input_header"],
                        "Sequence Length": results["sequence_length"],
                        "Class": get_class(matches[0][0]) if matches else "N/A",
                        "Best Genotype": matches[0][0] if matches else "N/A",
                        "Similarity %": matches[0][1] if matches else "N/A",
                        "Pathogenicity": cleavage["pathogenicity"]
                        if cleavage["cleavage_region_found"]
                        else "N/A",
                        "Motif": f"{cleavage['motif_type']} ({cleavage['motif_category']})"
                        if cleavage["cleavage_region_found"]
                        else "N/A",
                    }
                )

            export_df = pd.DataFrame(export_rows)
            st.dataframe(export_df, width="stretch", hide_index=True)


# ============================================================================
# TAB 2: PHYLOGENETIC TREES
# ============================================================================


# fonction très complexe pour rien qui m'a servi à personnalisé l'arbre comme ce que je voulais exactement
def tree_to_plotly(tree):

    # initialisation des listes / coordonnées
    x_lines, y_lines, x_nodes, y_nodes, labels = [], [], [], [], []
    y_pos = {}
    counter = [0]

    # fonction qui équilibre les espaces entre les feuilles de l'arbre
    def get_y(clade):
        if clade.is_terminal():
            y_pos[clade] = counter[0]
            counter[0] += 1
        else:
            for c in clade.clades:
                get_y(c)
            y_pos[clade] = sum(y_pos[c] for c in clade.clades) / len(clade.clades)

    # fonction qui assigne la profondeur dans l'arbre (racine = 0, parents = 1 ...)
    def get_x(clade, x=0):
        clade.x = x
        for c in clade.clades:
            get_x(c, x + 1)

    # calcul des positions
    get_y(tree.root)
    get_x(tree.root)

    # trouve la profondeur max de l'arbre
    max_x = max(c.x for c in tree.find_clades())

    # construit chaque trait de l'arbre pour l'affichage
    def collect(clade):
        for c in clade.clades:
            x_lines.append(clade.x)
            x_lines.append(c.x)
            x_lines.append(None)
            y_lines.append(y_pos[clade])
            y_lines.append(y_pos[clade])
            y_lines.append(None)
            x_lines.append(c.x)
            x_lines.append(c.x)
            x_lines.append(None)
            y_lines.append(y_pos[clade])
            y_lines.append(y_pos[c])
            y_lines.append(None)
            collect(c)

        # étent la ligne final pour afficher correctement le header
        if clade.is_terminal():
            x_lines.append(clade.x)
            x_lines.append(clade.x + 0.5)  # l'extension
            x_lines.append(None)
            y_lines.append(y_pos[clade])
            y_lines.append(y_pos[clade])
            y_lines.append(None)
            x_nodes.append(clade.x + 0.5)  # l'extension
            labels.append((clade.name or "")[:40])
        else:
            x_nodes.append(clade.x)  # position normal des nodes
            confidence = clade.confidence  # calcul et affiche les bootstraps
            if confidence is not None and confidence >= 0.5:
                labels.append(f"{confidence * 100:.1f}")
            else:
                labels.append("")

        y_nodes.append(y_pos[clade])

    # récupère toute les coordonnées pour ploty
    collect(tree.root)
    return x_lines, y_lines, x_nodes, y_nodes, labels, counter[0], max_x


with tab_tree:
    st.header("Phylogenetic Tree")

    if "all_sequences" not in st.session_state or "all_results" not in st.session_state:
        st.info("Run an analysis first to display the phylogenetic tree.")
    else:
        all_sequences = st.session_state["all_sequences"]

        # Construction de tous les arbres une seule fois et les mettre en cache
        if "trees" not in st.session_state:
            st.session_state["trees"] = {}

        for header, sequence in all_sequences.items():
            if header not in st.session_state["trees"]:
                with st.spinner(f"Building tree for {header[:40]}..."):
                    # trouve les 20 séquences les plus proche
                    neighbours = find_closest_neighbours(sequence, references, n=20)
                    # écrit un fasta temporaire pour les stocker
                    tmp_dir = os.path.join(DATA_FOLDER, "tmp")
                    os.makedirs(tmp_dir, exist_ok=True)
                    tmp_fasta = os.path.join(tmp_dir, "tmp_input.fasta")
                    tmp_aligned = os.path.join(tmp_dir, "tmp_aligned.fasta")
                    write_temp_fasta(header, sequence, neighbours, tmp_fasta)
                    # alignes ces séquences avec mafft
                    align_sequences_mafft(tmp_fasta, tmp_aligned)
                    # les sauvegardes dans session_state
                    st.session_state["trees"][header] = build_tree_fasttree(tmp_aligned)

        # Créer un onglet par séquence analysé
        tree_tabs = st.tabs([h[:30] for h in all_sequences.keys()])

        # L'affichage final
        for i, (header, sequence) in enumerate(all_sequences.items()):
            with tree_tabs[i]:
                tree = st.session_state["trees"][header]

                x_lines, y_lines, x_nodes, y_nodes, labels, n_leaves, max_x = (
                    tree_to_plotly(tree)
                )

                node_colors = []
                for label in labels:
                    if label.startswith("QUERY_"):
                        node_colors.append(
                            "#00FF00"
                        )  # la couleur de la séquence analysé
                    else:
                        node_colors.append(get_color(label))

                fig = go.Figure()
                # Les paramètres des lignes
                fig.add_trace(
                    go.Scatter(
                        x=x_lines,
                        y=y_lines,
                        mode="lines",
                        line=dict(color="rgba(150,150,150,0.3)", width=1.5),
                        hoverinfo="none",
                    )
                )
                # les paramètres du texte
                fig.add_trace(
                    go.Scatter(
                        x=x_nodes,
                        y=y_nodes,
                        mode="markers+text",
                        marker=dict(size=6, color=node_colors),
                        text=labels,
                        textposition="middle right",
                        textfont=dict(
                            size=12,
                            color=[
                                "#2ECC71" if lab.startswith("QUERY_") else "white"
                                for lab in labels
                            ],
                        ),
                        hoverinfo="text",
                    )
                )
                # les paramètres du graph
                fig.update_layout(
                    title=dict(
                        text=f"{header}",
                        font=dict(size=24, color="white"),
                        x=0,
                        xanchor="left",
                    ),
                    showlegend=False,
                    height=max(400, n_leaves * 25),
                    width=1800,
                    margin=dict(l=0, r=0, t=40, b=0),
                    xaxis=dict(
                        showgrid=False,
                        zeroline=False,
                        showticklabels=False,
                        fixedrange=False,
                    ),
                    yaxis=dict(
                        showgrid=False,
                        zeroline=False,
                        showticklabels=False,
                        fixedrange=False,
                    ),
                    plot_bgcolor="rgba(0,0,0,0)",
                    paper_bgcolor="rgba(0,0,0,0)",
                    shapes=[
                        dict(
                            type="rect",
                            xref="paper",
                            yref="paper",
                            x0=0,
                            y0=0,
                            x1=1,
                            y1=1,
                            line=dict(
                                color="rgba(255,255,255,0.3)", width=1
                            ),  # subtle white border
                        )
                    ],
                    dragmode="pan",
                )
                st.plotly_chart(
                    fig,
                    width="stretch",
                    config={"scrollZoom": True, "displayModeBar": True},
                )


# ============================================================================
# TAB 3: GLOBAL DISTRIBUTION MAP
# ============================================================================

with tab_map:
    st.header("Global Distribution Map")
    st.markdown("Explore the geographic distribution of NDV genotypes worldwide.")

    # chargement de toutes les données
    @st.cache_resource
    def load_locations_database():
        df_china = pd.read_csv(os.path.join(LOCATION_FOLDER, "china_provinces.csv"))
        df_us = pd.read_csv(os.path.join(LOCATION_FOLDER, "US_capitals.csv"))
        df_world = pd.read_csv(os.path.join(LOCATION_FOLDER, "countries.csv"))
        df_russia = pd.read_csv(os.path.join(LOCATION_FOLDER, "russia_regions.csv"))
        return df_china, df_us, df_world, df_russia

    df_china, df_us, df_world, df_russia = load_locations_database()

    # définition des centroïdes des gros pays
    CENTROIDS = {
        "USA": (37.09, -95.71),
        "China": (35.86, 104.19),
        "Russia": (61.52, 105.31),
    }

    # fonction qui retourne la lat, lon et les labels
    def parse_header_location(header):
        parts = header.split("_")
        padded = "_" + "_".join(parts) + "_"  # pour chercher _XX_ proprement

        # Détection du pays
        country = None
        if "_USA_" in padded:
            country = "USA"
        elif "_China_" in padded:
            country = "China"
        elif "_Russia_" in padded:
            country = "Russia"
        else:
            # Cherche dans la colonne country du CSV monde
            for _, row in df_world.iterrows():
                if f"_{row['country']}_" in padded:
                    return row["lat"], row["lon"], row["country"]
            return None, None, "Unknown"

        # Détection de la région
        if country == "USA":
            # abréviation (_PA_, _NY_...)
            for _, row in df_us.iterrows():
                if f"_{row['abbreviation_us']}_" in padded:
                    return row["lat"], row["lon"], f"{row['state']}, USA"
            # nom complet (_New_York_...)
            for _, row in df_us.iterrows():
                state_fmt = row["state"]  # déjà avec _ dans le CSV ex: New_York
                if f"_{state_fmt}_" in padded:
                    return row["lat"], row["lon"], f"{row['state']}, USA"
            # Fallback
            return CENTROIDS["USA"][0], CENTROIDS["USA"][1], "USA"

        if country == "China":
            # abréviation (_AH_, _BJ_...)
            for _, row in df_china.iterrows():
                if f"_{row['abbreviation_cn']}_" in padded:
                    return row["lat"], row["lon"], f"{row['province']}, China"
            # nom complet (_Shanghai_...)
            for _, row in df_china.iterrows():
                if f"_{row['province']}_" in padded:
                    return row["lat"], row["lon"], f"{row['province']}, China"
            # Fallback
            return CENTROIDS["China"][0], CENTROIDS["China"][1], "China"

        if country == "Russia":
            # Nom de région (_Novosibirsk_, _FarEast_, _Amur_region_...)
            for _, row in df_russia.iterrows():
                if f"_{row['region_ru']}_" in padded:
                    return row["lat"], row["lon"], f"{row['region_ru']}, Russia"
            # Fallback
            return CENTROIDS["Russia"][0], CENTROIDS["Russia"][1], "Russia"

    SEQUENCES_FOLDER = os.path.join(DATA_FOLDER, "sequences")

    # création d'un dataframe avec header, loc et genotype
    @st.cache_resource
    def build_map_dataframe():
        rows = []
        for file in os.listdir(SEQUENCES_FOLDER):
            if file.endswith(".fas"):
                path = os.path.join(SEQUENCES_FOLDER, file)
                with open(path, "r") as f:
                    for line in f:
                        if line.startswith(">"):
                            header = line.strip().lstrip(">")
                            parts = header.split("_")  # récuperation des génotypes
                            genotype = parts[1] if len(parts) > 1 else "Unknown"
                            lat, lon, label = parse_header_location(header)
                            rows.append(
                                {
                                    "header": header,
                                    "genotype": genotype,
                                    "lat": lat,
                                    "lon": lon,
                                    "label": label,
                                }
                            )
        return pd.DataFrame(rows)

    df_map = build_map_dataframe()

    # Nettoyer le DataFrame
    df_map_clean = df_map[
        (df_map["label"] != "Unknown")
        & (df_map["lat"].notna())
        & (df_map["lon"].notna())
    ].copy()

    # Déterminer la classe de chaque génotype
    def get_class(genotype):
        if genotype.startswith("I.") or genotype == "I":
            return "Class I"
        else:
            return "Class II"

    df_map_clean["class"] = df_map_clean["genotype"].apply(get_class)

    # Affichage des filtres + Map
    col_filters, col_map = st.columns([1, 5])

    available_genotypes = sorted(df_map_clean["genotype"].unique())

    with col_filters:
        st.subheader("Filters")
        selected_genotypes = st.multiselect(
            "Genotypes",
            options=available_genotypes,
            default=[],
        )

    # Appliquer les filtres
    df_filtered = df_map_clean[df_map_clean["genotype"].isin(selected_genotypes)]

    df_agg = (
        df_filtered.groupby(["lat", "lon", "label"]).size().reset_index(name="count")
    )

    # affichage de la carte avec pydeck
    with col_map:
        st.pydeck_chart(
            pdk.Deck(
                map_style="https://basemaps.cartocdn.com/gl/dark-matter-gl-style/style.json",
                initial_view_state=pdk.ViewState(
                    latitude=30.0,
                    longitude=20.0,
                    zoom=1,
                    pitch=20,
                ),
                layers=[
                    pdk.Layer(
                        "HeatmapLayer",
                        data=df_filtered,
                        get_position="[lon, lat]",
                        opacity=0.6,
                        radius_pixels=70,
                        color_range=[
                            [0, 99, 153, 120],  # rappel le 4ème nombres est l'intensité
                            [0, 120, 170, 180],
                            [0, 153, 180, 190],
                            [0, 201, 167, 210],
                            [100, 220, 180, 240],
                            [255, 255, 200, 255],
                        ],
                    ),
                    pdk.Layer(
                        "ScatterplotLayer",
                        data=df_agg,
                        get_position="[lon, lat]",
                        get_fill_color=[0, 0, 0, 200],
                        get_line_color=[100, 220, 180, 255],
                        stroked=True,
                        line_width_min_pixels=1,
                        get_radius=30000,
                        pickable=True,
                    ),
                ],
                tooltip={"text": "{label}\nSequences: {count}"},
            )
        )


# ============================================================================
# TAB 4: Help
# ============================================================================

with help_tab:
    st.divider()
    st.markdown(
        "##### **If you encounter any issues feel free to report them [here](https://github.com/Soupirr/NDV-genotyper/issues).**"
    )
    st.divider()

    info_tab, stat_tab = st.tabs(["Informations", "Statistics"])

    with info_tab:
        with open("QUICK_START.md", "r", encoding="utf-8") as f:
            informations = f.read()

        st.markdown(informations)

    with stat_tab:
        st.header("Database Statistics")

        db_references, db_files_count, db_total_count, _ = load_all_references()

        # analyse de tout les champs de chaque header
        genotype_counts = {}
        class_counts = {"Class I": 0, "Class II": 0}
        host_counts = {}
        country_counts = {}
        year_counts = {}

        # Normalisation des noms d'hôtes
        # Pour ajouter un hôte rajouter une ligne raw_host,normalized_host  dans le CSV
        _df_host_map = pd.read_csv(
            os.path.join(DATA_FOLDER, "hosts", "host_normalize.csv")
        )
        HOST_NORMALIZE = dict(
            zip(_df_host_map["raw_host"].str.lower(), _df_host_map["normalized_host"])
        )

        KNOWN_COUNTRIES = set(df_world["country"].str.lower()) | {
            "china",
            "usa",
            "russia",
        }

        # US abbreviations
        US_STATE_ABBREVS = set(df_us["abbreviation_us"].str.lower())

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

            # Filtrage des pays
            country_idx = host_idx + 1
            if len(parts) > country_idx:
                raw_country = parts[country_idx]
                if raw_country and not raw_country.isdigit() and len(raw_country) > 1:
                    country_counts[raw_country] = country_counts.get(raw_country, 0) + 1

            # Filtrage des années
            last = parts[-1]
            if last.isdigit() and len(last) == 4:
                year_counts[last] = year_counts.get(last, 0) + 1

        df_genotypes = pd.DataFrame(
            sorted(genotype_counts.items(), key=lambda x: x[1], reverse=True),
            columns=["Genotype", "Sequences"],
        )

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

            # palette importé depuis ColorBrewer
            PALETTE = [
                "#0099cc",
                "#00c9a7",
                "#f0a500",
                "#e05c5c",
                "#7b68ee",
                "#20b2aa",
                "#ff7f50",
                "#9acd32",
                "#ba55d3",
                "#4682b4",
                "#cd853f",
                "#5f9ea0",
                "#d2691e",
                "#6495ed",
                "#dc143c",
            ]
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
