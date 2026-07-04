"""Phylogenetic Tree tab."""

import os
import time
import streamlit as st
import plotly.graph_objects as go
from genotyper.analyzer import (
    find_closest_neighbours,
    write_temp_fasta,
    align_sequences_mafft,
    build_tree_fasttree,
    build_tree_iqtree2,
    tree_to_newick,
    get_color,
    clean_sequence,
)
from genotyper.config import DATA_FOLDER, PALETTE
from genotyper import report
from genotyper.tabs.analyze_tab import load_all_references

# ============================================================================


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
            if st.session_state.get("tree_mode") == "Phylogram":
                get_x(c, x + (c.branch_length or 0))
            else:
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
            labels.append((clade.name or "")[:70])
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


def draw_tree(tree, title, multi_query=False):
    x_lines, y_lines, x_nodes, y_nodes, labels, n_leaves, max_x = tree_to_plotly(tree)

    if multi_query:
        q_labels = list(
            dict.fromkeys(lab for lab in labels if lab.startswith("QUERY_"))
        )
        q_colors = {lab: PALETTE[i % len(PALETTE)] for i, lab in enumerate(q_labels)}
    else:
        q_colors = {}

    node_colors = [
        q_colors.get(lab, "#00FF00") if lab.startswith("QUERY_") else get_color(lab)
        for lab in labels
    ]

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
                    q_colors.get(lab, "#2ECC71")
                    if lab.startswith("QUERY_")
                    else "white"
                    for lab in labels
                ],
            ),
            hoverinfo="text",
        )
    )
    # les paramètres du graph
    fig.update_layout(
        title=dict(
            text=f"{title[:66]}",
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
        plot_bgcolor="#060d14",
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
    # stash the figure so the report can embed it (identical to what's shown)
    st.session_state.setdefault("report_tree_figs", {})[title] = fig
    st.plotly_chart(
        fig,
        width="stretch",
        config={"scrollZoom": True, "displayModeBar": True},
    )
    st.download_button(
        "Download tree (Newick)",
        data=tree_to_newick(tree),
        file_name=f"{title[:50].replace(' ', '_')}.nwk",
        mime="text/plain",
        key=f"newick_dl_{title}",
    )


# ============================================================================


def render(path):
    references, _, _, _ = load_all_references(path)
    st.header("Phylogenetic Tree")

    if "all_sequences" not in st.session_state or "all_results" not in st.session_state:
        st.info("Run an analysis first to display the phylogenetic tree.")
    elif "load_tree" not in st.session_state:
        col_slider_tree, col_button_tree = st.columns([3, 1])
        with col_slider_tree:
            tree_type = st.radio(
                "Tree type",
                ["Combined", "Per-query"],
                horizontal=True,
                help="Per-query: one tree each. Combined: all queries in a single tree.",
            )
            if tree_type == "Combined":
                n_per_query = st.slider(
                    "Neighbours per query",
                    0,
                    20,
                    5,
                    help="Closest references added per query (keep low to stay readable).",
                )
            else:
                n_neighbours = st.slider("Number of Neighbours", 5, 100, 20, step=5)
            tree_mode = st.radio(
                "Tree Mode",
                ["Cladogram", "Phylogram"],
                horizontal=True,
                help="Cladogram: uniform branch lengths. Phylogram: real evolutionary distances.",
            )

            tree_method = st.radio(
                "Tree Method",
                ["FastTree", "IQ-TREE2"],
                horizontal=True,
                help="FastTree: fast, approximate. IQ-TREE2: full ML + ModelFinder + bootstrap, slower for publication.",
            )
            if tree_method == "IQ-TREE2":
                st.caption(
                    "IQ-TREE2 is significantly slower than FastTree, especially in Per-query mode with a large number of sequences."
                )
        with col_button_tree:
            st.write("")  # petit espace pour aligner verticalement avec le slider
            if st.button("Create Trees", type="primary"):
                st.session_state["load_tree"] = True
                st.session_state["tree_type"] = tree_type
                st.session_state["tree_mode"] = tree_mode
                st.session_state["tree_method"] = tree_method
                if tree_type == "Combined":
                    st.session_state["n_per_query"] = n_per_query
                else:
                    st.session_state["n_neighbours"] = n_neighbours
                st.rerun()

    else:
        all_sequences = st.session_state["all_sequences"]
        tree_type = st.session_state.get("tree_type", "Per-query")
        tmp_dir = os.path.join(DATA_FOLDER, "tmp")
        os.makedirs(tmp_dir, exist_ok=True)

        if tree_type == "Combined":
            if "combined_tree" not in st.session_state:
                with st.spinner("Building combined tree..."):
                    n_per = st.session_state["n_per_query"]
                    pool = {}
                    for qseq in all_sequences.values():
                        for h, s, _ in find_closest_neighbours(
                            qseq, references, n=n_per
                        ):
                            pool[h] = s

                    tmp_fasta = os.path.join(tmp_dir, "combined_input.fasta")
                    tmp_aligned = os.path.join(tmp_dir, "combined_aligned.fasta")
                    with open(tmp_fasta, "w") as f:
                        for qh, qs in all_sequences.items():
                            f.write(f">QUERY_{qh}\n{clean_sequence(qs)}\n")
                        for h, s in pool.items():
                            f.write(f">{h}\n{clean_sequence(s)}\n")

                    align_sequences_mafft(tmp_fasta, tmp_aligned)
                    tree_builder = (
                        build_tree_iqtree2
                        if st.session_state.get("tree_method") == "IQ-TREE2"
                        else build_tree_fasttree
                    )
                    st.session_state["combined_tree"] = tree_builder(tmp_aligned)

            tree = st.session_state["combined_tree"]
            if tree is None:
                st.warning(
                    "Could not build the combined tree — try re-running the analysis."
                )
            else:
                draw_tree(tree, "Combined tree — all queries", multi_query=True)

        # Construction de tous les arbres une seule fois et les mettre en cache
        else:
            if "trees" not in st.session_state:
                st.session_state["trees"] = {}

            for header, sequence in all_sequences.items():
                if header not in st.session_state["trees"]:
                    with st.spinner(f"Building tree for {header[:40]}..."):
                        # trouve les 20 séquences les plus proche
                        neighbours = find_closest_neighbours(
                            sequence, references, n=st.session_state["n_neighbours"]
                        )
                        # écrit un fasta temporaire pour les stocker
                        tmp_fasta = os.path.join(tmp_dir, "tmp_input.fasta")
                        tmp_aligned = os.path.join(tmp_dir, "tmp_aligned.fasta")
                        write_temp_fasta(header, sequence, neighbours, tmp_fasta)
                        # alignes ces séquences avec mafft
                        align_sequences_mafft(tmp_fasta, tmp_aligned)
                        # les sauvegardes dans session_state
                        tree_builder = (
                            build_tree_iqtree2
                            if st.session_state.get("tree_method") == "IQ-TREE2"
                            else build_tree_fasttree
                        )
                        tree = tree_builder(tmp_aligned)
                        if tree is None:
                            st.warning(
                                f"Could not build tree for {header[:40]} — try re-running the analysis."
                            )
                        else:
                            st.session_state["trees"][header] = tree

            # Créer un onglet par séquence analysé
            tree_tabs = st.tabs([h[:30] for h in all_sequences.keys()])

            # L'affichage final
            for i, (header, sequence) in enumerate(all_sequences.items()):
                with tree_tabs[i]:
                    if header not in st.session_state["trees"]:
                        st.warning(f"Tree unavailable for {header[:40]}")
                        continue
                    tree = st.session_state["trees"][header]
                    draw_tree(tree, header)

    # --- Export Report ---
    if "report_export_df" in st.session_state:
        st.divider()
        st.subheader("Export Report")
        tree_figs = st.session_state.get("report_tree_figs", {})
        matrix_fig = st.session_state.get("report_matrix_fig")
        all_results = st.session_state.get("all_results", [])
        method = st.session_state.get("method", "hamming")
        elapsed_time = st.session_state.get("elapsed_time", 0)
        export_df = st.session_state["report_export_df"]

        if not tree_figs:
            st.caption("Tip: build trees above first to include them in the report.")
        if st.button("Generate Report", type="primary"):
            st.session_state["report_html"] = report.build_report_html(
                all_results=all_results,
                export_df=export_df,
                method=method,
                elapsed_time=elapsed_time,
                matrix_fig=matrix_fig,
                tree_figs=tree_figs,
                entry_name=os.path.basename(path),
            )
        if "report_html" in st.session_state:
            st.download_button(
                "Download Report (HTML)",
                data=st.session_state["report_html"],
                file_name=f"genotyper_report_{time.strftime('%Y%m%d_%H%M%S', time.localtime())}.html",
                mime="text/html",
            )
