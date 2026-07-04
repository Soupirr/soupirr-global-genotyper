"""Analyze Sequences tab."""

import os
import time
import platform
import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from genotyper.analyzer import (
    FASTAParser,
    analyze_sequence,
    unpack_top_match,
    SequenceSimilarity,
)


if platform.system() == "Windows":
    import winsound

# ============================================================================

# Charge les datasets de référence / reset les dictionnaires
_FASTA_EXTS = {".fasta", ".fas", ".fa", ".txt"}


# Charge les datasets de référence / reset les dictionnaires
@st.cache_resource
def load_all_references(path):
    combined_sequences = {}
    errors = []

    fas_files = [
        f for f in os.listdir(path) if os.path.splitext(f)[1].lower() in _FASTA_EXTS
    ]

    for filename in fas_files:
        try:
            seqs = FASTAParser.parse_file(os.path.join(path, filename))
            if seqs:
                combined_sequences.update(seqs)
        except Exception as e:
            errors.append(f"{filename}: {e}")

    return combined_sequences, len(fas_files), len(combined_sequences), errors


# ============================================================================


def render(path, config=None):
    references, files_count, total_count, load_errors = load_all_references(path)
    pathogenicity_config = config or {}

    st.header(
        f"Sequence Analyzer Tool \u00a0\u00a0\u00a0({total_count} Sequences Loaded)",
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
        1. Download the reference datasets from: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/
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
            ["Pairwise", "Hamming"],
            help="Hamming: simple mismatch count. Pairwise: accounts for insertions/deletions",
        )

        # Alerte pour la méthode de Leveshtein
        if "Hamming" in similarity_method:
            st.warning(
                "Faster but only works if sequences are only aligned to references."
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
            st.session_state.pop("report_matrix_fig", None)
            st.session_state.pop("report_export_df", None)

        with st.form("analyze_form", border=False):
            input_tab1, input_tab2 = st.tabs(["Paste FASTA", "Upload File"])

            # Input sous forme de texte
            with input_tab1:
                input_fasta = st.text_area(
                    "Paste your sequence in FASTA format:",
                    height=200,
                    placeholder=">your_sequence_name\nATGGGCTCCAGATCCTCTAC...",
                    help="Should be a complete gene sequence",
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
            st.session_state.pop("report_matrix_fig", None)
            st.session_state.pop("report_export_df", None)
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
                    result = analyze_sequence(
                        input_fasta=single_fasta,
                        reference_sequences=references,
                        top_matches=top_n_matches,
                        similarity_method=method,
                        pathogenicity_config=pathogenicity_config,
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

        st.divider()

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
                    col1, col2, col3, col4, col5 = st.columns(5)
                    with col1:
                        st.metric("Best Match Genotype", top_match["genotype"])
                    with col2:
                        st.metric("Reference Sequences", top_match["sample_count"])
                    with col3:
                        st.metric("Avg Match Score", f"{top_match['avg_similarity']}%")
                    with col4:
                        st.metric("Best Match Score", f"{top_match['best_score']}%")
                    with col5:
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
                            "⚠ Average Similarity score is below 95% - result may not be interpretable with the Hamming method. "
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
                            "Best Match Header",
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

                if not results.get("pathogenicity_configured"):
                    st.info(
                        "No pathogenicity configuration for this entry. "
                        "Add a cleavage site position and motif file when creating the entry to enable this analysis."
                    )
                else:

                    def display_reading_frame(cleavage: dict):
                        if cleavage["cleavage_region_found"]:
                            motif = cleavage["motif_type"]
                            protein = cleavage["cleavage_protein"]
                            motif_found = bool(motif and protein and motif in protein)

                            col1, col2, col3 = st.columns([3, 1, 1])

                            with col1:
                                st.write("**Cleavage Region (protein):**")
                                if motif_found:
                                    idx = protein.index(motif)
                                    before = protein[:idx]
                                    after = protein[idx + len(motif) :]
                                    color = color = "#00dac7"
                                    st.markdown(
                                        f"<div style='font-size:20px; text-align:center; font-family:monospace; background-color:#0e1117; padding:8px; border-radius:4px; border:1px solid rgba(255,255,255,0.2);'>{before}<span style='color:{color}'>{motif}</span>{after}</div>",
                                        unsafe_allow_html=True,
                                    )
                                else:
                                    st.markdown(
                                        f"<div style='font-size:20px; text-align:center; font-family:monospace; background-color:#0e1117; padding:8px; border-radius:4px; border:1px solid rgba(255,255,255,0.2); color:#868d2f;'>{protein}</div>",
                                        unsafe_allow_html=True,
                                    )
                            with col2:
                                st.write("**Motif Found:**")
                                if motif_found:
                                    st.code(motif, language="text")
                                else:
                                    st.code("Not found", language="text")
                            with col3:
                                st.write("**Motif Category:**")
                                st.code(
                                    f"{cleavage['motif_category'] or 'Unknown'}",
                                    language="text",
                                )

                            pathogenicity = cleavage["pathogenicity"]

                            if pathogenicity == "Undetermined":
                                st.warning("Undetermined")
                                st.write(
                                    "No known cleavage site motif was found in the expected region. "
                                    "The sequence may be incomplete or represent an unusual strain."
                                )
                            else:
                                st.info(f"**{pathogenicity.capitalize()}**")
                        else:
                            st.warning(
                                "Could not analyze cleavage site. Sequence may be incomplete."
                            )

                    cadre_1, cadre_2, cadre_3 = st.tabs(
                        ["Main Reading Frame", "Reading Frame +1", "Reading Frame -1"]
                    )
                    with cadre_1:
                        display_reading_frame(results["cleavage_main"])
                    with cadre_2:
                        display_reading_frame(results["cleavage_plus_one"])
                    with cadre_3:
                        display_reading_frame(results["cleavage_minus_one"])

                    st.caption(
                        "Pathogenicity prediction is based on cleavage site motif matching and should be "
                        "interpreted as indicative only. Results must be validated by biological assays "
                        "before any clinical or regulatory conclusion."
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
            st.session_state["report_matrix_fig"] = fig_matrix

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

            top = unpack_top_match(matches[0]) if matches else None
            export_rows.append(
                {
                    "Header": results["input_header"],
                    "Length": results["sequence_length"],
                    "Best Genotype": top["genotype"] if top else "N/A",
                    "Ref Seq": top["sample_count"] if top else "N/A",
                    "Avg Match Score (%)": top["avg_similarity"] if top else "N/A",
                    "Best Match Score (%)": top["best_score"] if top else "N/A",
                    "Pathogenicity": best_cleavage["pathogenicity"]
                    if best_cleavage
                    else "Undetermined",
                    "Cleavage Motif": best_cleavage["motif_type"]
                    if best_cleavage
                    else "N/A",
                    "Motif Category": best_cleavage["motif_category"]
                    if best_cleavage
                    else "N/A",
                    "Reading Frame": best_frame_label,
                }
            )

        export_df = pd.DataFrame(export_rows)
        st.dataframe(export_df, width="stretch", hide_index=True)
        st.session_state["report_export_df"] = export_df
