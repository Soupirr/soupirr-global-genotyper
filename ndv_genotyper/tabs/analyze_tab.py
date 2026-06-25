"""Analyze Sequences tab."""

import os
import time
import platform
import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from ndv_genotyper.analyzer import (
    FASTAParser,
    analyze_newcastle_sequence,
    unpack_top_match,
    get_class,
    SequenceSimilarity,
)
from ndv_genotyper.config import SEQ_FOLDER
from ndv_genotyper import report

if platform.system() == "Windows":
    import winsound

# ============================================================================


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

# ============================================================================


def render():
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
                "Analysis can take **5-30+ minutes** depending on number of sequences and computer power."
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
            st.session_state.pop("combined_tree", None)
            st.session_state.pop("report_tree_figs", None)
            st.session_state.pop("report_html", None)

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
            st.session_state.pop("combined_tree", None)
            st.session_state.pop("report_tree_figs", None)
            st.session_state.pop("report_html", None)
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
        st.write("")
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
                            "Hamming" if "Hamming" in similarity_method else "Pairwise",
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
                            st.success(f"Low-virulence (Confidence: {confidence_path})")
                            st.write(
                                "This sequence shows characteristics of a low-virulence strain with monobasic or L117 cleavage site."
                            )
                        else:
                            st.warning(f"Undetermined (Confidence: {confidence_path})")
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
                height=700,
                width=850,
                plot_bgcolor="#060d14",
                paper_bgcolor="#060d14",
                margin=dict(t=40, l=50, r=50, b=50),
            )
            st.plotly_chart(fig_matrix, width="content")

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
                    "Reading Frame": best_frame_label if best_frame_label else "N/A",
                }
            )

        export_df = pd.DataFrame(export_rows)
        st.dataframe(export_df, width="stretch", hide_index=True)

        # --- Full HTML report ---
        st.divider()
        st.subheader("Export Report")
        tree_figs = st.session_state.get("report_tree_figs", {})
        matrix_fig = fig_matrix if (len(all_sequences) > 1 and show_matrix) else None
        if not tree_figs:
            st.caption(
                "Tip: build trees in the Phylogenetic Tree tab first to include them in the report."
            )
        if st.button("Generate Report", type="primary"):
            st.session_state["report_html"] = report.build_report_html(
                all_results=all_results,
                export_df=export_df,
                method=method,
                elapsed_time=elapsed_time,
                matrix_fig=matrix_fig,
                tree_figs=tree_figs,
            )
        if "report_html" in st.session_state:
            st.download_button(
                "Download Report (HTML)",
                data=st.session_state["report_html"],
                file_name=f"NDV_report_{time.strftime('%Y%m%d_%H%M%S', time.localtime())}.html",
                mime="text/html",
            )
