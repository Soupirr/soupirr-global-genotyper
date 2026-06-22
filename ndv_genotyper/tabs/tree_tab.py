"""Phylogenetic Tree tab."""

import os
import streamlit as st
import plotly.graph_objects as go
from ndv_genotyper.analyzer import (
    find_closest_neighbours,
    write_temp_fasta,
    align_sequences_mafft,
    build_tree_fasttree,
    get_color,
)

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_FOLDER = os.path.join(_ROOT, "data")

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------


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


# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------


def render(references):
    st.header("Phylogenetic Tree")

    if "all_sequences" not in st.session_state or "all_results" not in st.session_state:
        st.info("Run an analysis first to display the phylogenetic tree.")
    elif "load_tree" not in st.session_state:
        col_slider_tree, col_button_tree = st.columns([3, 1])
        with col_slider_tree:
            n_neighbours = st.slider(
                "Number of Neighbours",
                min_value=5,
                max_value=100,
                value=20,
                step=5,
                help="Number of closest reference sequences used to build the phylogenetic tree (default: 20)",
            )
            tree_mode = st.radio(
                "Tree Mode",
                ["Cladogram", "Phylogram"],
                horizontal=True,
                help="Cladogram: uniform branch lengths. Phylogram: real evolutionary distances.",
            )
        with col_button_tree:
            st.write("")  # petit espace pour aligner verticalement avec le slider
            if st.button("Create Trees", type="primary"):
                st.session_state["load_tree"] = True
                st.session_state["n_neighbours"] = n_neighbours
                st.session_state["tree_mode"] = tree_mode
                st.rerun()

    else:
        all_sequences = st.session_state["all_sequences"]

        # Construction de tous les arbres une seule fois et les mettre en cache
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
                    tmp_dir = os.path.join(DATA_FOLDER, "tmp")
                    os.makedirs(tmp_dir, exist_ok=True)
                    tmp_fasta = os.path.join(tmp_dir, "tmp_input.fasta")
                    tmp_aligned = os.path.join(tmp_dir, "tmp_aligned.fasta")
                    write_temp_fasta(header, sequence, neighbours, tmp_fasta)
                    # alignes ces séquences avec mafft
                    align_sequences_mafft(tmp_fasta, tmp_aligned)
                    # les sauvegardes dans session_state
                    tree = build_tree_fasttree(tmp_aligned)
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
                        text=f"{header[:66]}",
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
