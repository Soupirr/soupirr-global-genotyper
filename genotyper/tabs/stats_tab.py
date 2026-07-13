"""Statistics / Help tab."""

import os
import re
import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from genotyper.tabs import map_tab, analyze_tab
from genotyper.analyzer import CleavageSiteAnalyzer
from genotyper.config import PALETTE, PATHO_PALETTE

# ============================================================================


def render(path, entry_config=None):
    st.write("")

    if "data_loaded" not in st.session_state:
        st.info(
            "Building the statistics scans every sequence in the database. "
            "Click to load."
        )
        if st.button("Load Statistics", type="primary"):
            st.session_state["data_loaded"] = True
            st.rerun()
        return

    entry_name = os.path.basename(path)
    md_files = [f for f in os.listdir(path) if f.endswith(".md")]
    has_patho_md = len(md_files) > 0
    patho_md_path = os.path.join(path, md_files[0]) if has_patho_md else None

    tab_labels = ["Statistics", "Map"]
    if has_patho_md:
        tab_labels.append("Pathogenicity Documentation")

        # Détection multi/mono avant l'ouverture des tabs
    db_references, db_files_count, db_total_count, _ = analyze_tab.load_all_references(
        path
    )
    is_multi = bool(entry_config and entry_config.get("multi", False))

    # Métriques globales avant le graphe et la sélection de gène
    if is_multi:
        _first_refs = db_references[list(db_references.keys())[0]]
        _total_sequences = sum(len(refs) for refs in db_references.values())
        _total_bp = sum(
            len(seq)
            for gene_seqs in db_references.values()
            for seq in gene_seqs.values()
        )
        _unique_genotypes = len(
            {
                h.split("|")[2]
                for gene_refs in db_references.values()
                for h in gene_refs
                if len(h.split("|")) >= 3
            }
        )
    else:
        _total_sequences = len(db_references)
        _total_bp = sum(len(seq) for seq in db_references.values())
        _unique_genotypes = len(
            {h.split("|")[2] for h in db_references if len(h.split("|")) >= 3}
        )

    if is_multi:
        col_m1, col_m2, col_m3 = st.columns(3)
        col_m1.metric("Total Sequences", f"{_total_sequences:,}")
        col_m2.metric("Total Base Pairs", f"{_total_bp:,}")
        col_m3.metric("Unique Genotypes", f"{_unique_genotypes:,}")
        st.divider()

    if is_multi:
        # Graphe global AVANT la sélection du gène
        gene_geno_counts = {gene: {} for gene in db_references}
        global_geno_totals = {}
        for gene, gene_refs in db_references.items():
            for header in gene_refs.keys():
                parts = header.split("|")
                if len(parts) < 7:
                    continue
                geno = parts[2]
                gene_geno_counts[gene][geno] = gene_geno_counts[gene].get(geno, 0) + 1
                global_geno_totals[geno] = global_geno_totals.get(geno, 0) + 1

        sorted_genos = [
            g
            for g, _ in sorted(
                global_geno_totals.items(), key=lambda x: x[1], reverse=True
            )
        ]

        st.subheader("All Genotypes (combined)")
        fig_global = go.Figure()
        gene_colors = ["#0099cc", "#00c9a7", "#f4a261", "#e76f51", "#8ecae6"]
        for i, gene in enumerate(db_references.keys()):
            counts = gene_geno_counts[gene]
            y_vals = [counts.get(g, 0) for g in sorted_genos]
            fig_global.add_trace(
                go.Bar(
                    name=gene,
                    x=sorted_genos,
                    y=y_vals,
                    marker=dict(
                        color=gene_colors[i % len(gene_colors)],
                    ),
                    text=y_vals,
                    textposition="inside",
                    hovertemplate="%{x}<br>%{y} seq<br>" + gene + "<extra></extra>",
                )
            )
        fig_global.update_layout(
            barmode="stack",
            xaxis_title="Genotype",
            yaxis_title="Nb. of Sequences",
            height=350,
            margin=dict(l=0, r=0, t=0, b=0),
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            font=dict(color="white"),
            xaxis=dict(type="category", tickangle=-45),
            legend=dict(orientation="h", y=1.1, x=0),
        )
        st.plotly_chart(fig_global, width="stretch")
        st.divider()

        # Sélection du gène APRÈS le graphe global
        gene_names = list(db_references.keys())
        selected_gene = st.segmented_control("Gene", gene_names, default=gene_names[0])
        if selected_gene not in gene_names:
            selected_gene = gene_names[0]
        flat_refs = db_references[selected_gene]
        gene_path = os.path.join(path, selected_gene)
        gene_cfg = (entry_config.get("genes", {}) or {}).get(selected_gene, {})
        db_total_count = len(flat_refs)
    else:
        flat_refs = db_references
        gene_path = path
        gene_cfg = entry_config or {}

    tabs = st.tabs(tab_labels)
    stat_tab, map_subtab = tabs[0], tabs[1]

    if has_patho_md:
        with tabs[2]:
            with open(patho_md_path, "r", encoding="utf-8") as p:
                st.markdown(p.read())

    with map_subtab:
        map_tab.render(gene_path)  # ← gene_path disponible maintenant

    with stat_tab:
        st.header("Database Statistics")

        # Métriques générales
        total_bp = sum(len(seq) for seq in flat_refs.values())
        unique_genotypes = len(
            {h.split("|")[2] for h in flat_refs if len(h.split("|")) >= 3}
        )

        col_m1, col_m2, col_m3 = st.columns(3)
        col_m1.metric("Total Sequences", f"{db_total_count:,}")
        col_m2.metric("Total Base Pairs", f"{total_bp:,}bp")
        col_m3.metric("Average Sequence Lenght", f"{(total_bp / db_total_count):.1f}bp")
        st.divider()

        # ── Comptage des sous-génotypes pour le gène sélectionné ──
        genotype_counts = {}
        host_counts = {}
        year_counts = {}

        # Pattern du gène sélectionné pour extraire le sous-génotype
        geno_pattern = gene_cfg.get("genotype_pattern", "") if is_multi else ""

        for header in flat_refs.keys():
            parts = header.split("|")
            if len(parts) < 7:
                continue

            raw_geno = parts[2]
            if geno_pattern:
                m = re.search(geno_pattern, raw_geno)
                geno = m.group(0) if m else raw_geno
            else:
                geno = raw_geno
            genotype_counts[geno] = genotype_counts.get(geno, 0) + 1

            host = parts[3].replace("_", " ")
            host = host if host not in ("UNKNOWN", "?", "") else "Unspecified"
            host_counts[host] = host_counts.get(host, 0) + 1

            year = parts[6]
            if year.isdigit() and len(year) == 4:
                year_counts[year] = year_counts.get(year, 0) + 1

        # Pour ajouter un hôte rajouter une ligne raw_host,normalized_host  dans le CSV
        patho_csv_path = os.path.join(path, f"{entry_name}_pathogenicity.csv")
        df_patho = (
            pd.read_csv(patho_csv_path) if os.path.exists(patho_csv_path) else None
        )

        df_genotypes = pd.DataFrame(
            sorted(genotype_counts.items(), key=lambda x: x[1], reverse=True),
            columns=["Genotype", "Sequences"],
        )

        _df_map_stats = map_tab.build_map_dataframe(gene_path)
        _df_map_stats = _df_map_stats[_df_map_stats["label"] != "Unknown"].copy()
        _df_map_stats["country"] = _df_map_stats["label"].apply(
            lambda x: x.split(",")[-1].strip()
        )
        country_counts = _df_map_stats["country"].value_counts().to_dict()

        # Statistique par gène sélectionné
        years_with_data = [int(y) for y in year_counts if year_counts[y] > 0]
        col1, col2, col3 = st.columns(3)
        col2.metric("Unique Sub-Genotypes", len(genotype_counts))
        col3.metric(
            "Year Range",
            f"{min(years_with_data)}–{max(years_with_data)}"
            if years_with_data
            else "N/A",
        )
        col1.metric("Unique Genotypes", unique_genotypes)

        st.divider()

        sub_label = (
            f"Sequences per Sub-Genotype — {selected_gene}"
            if is_multi
            else "Sequences per Genotype"
        )
        st.subheader(sub_label)
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
            xaxis=dict(type="category", tickangle=-45),
        )
        st.plotly_chart(fig_bar, width="stretch")

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
        df_health["% of Total"] = (df_health["Sequences"] / db_total_count * 100).map(
            "{:.1f}%".format
        )
        df_health = df_health[["Genotype", "Sequences", "% of Total", "Status"]]

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

        # Fonction qui permet de faire une analyse entière de la pathogénéicité du dataset de ref
        if df_patho is None:
            cfg = gene_cfg
            if "cleavage_start" not in cfg:
                st.info("No pathogenicity configuration for this entry.")
            else:
                st.info("No pathogenicity data generated yet for this entry.")
                if st.button("Generate Pathogenicity Data", type="primary"):
                    with st.spinner("Analyzing all reference sequences..."):
                        rows = []
                        for header, sequence in flat_refs.items():
                            parts = header.split("|")
                            genotype = parts[2] if len(parts) > 2 else "Unknown"
                            cleavage, _, _ = CleavageSiteAnalyzer.analyze(
                                sequence,
                                cleavage_start=cfg["cleavage_start"],
                                motifs_by_type=cfg.get("motifs_by_type"),
                            )
                            rows.append(
                                {
                                    "Header": header,
                                    "Best Genotype": genotype,
                                    "Pathogenicity": cleavage["pathogenicity"],
                                    "Cleavage Motif": cleavage["motif_type"] or "N/A",
                                    "Motif Category": cleavage["motif_category"]
                                    or "N/A",
                                }
                            )
                        df_generated = pd.DataFrame(rows)
                        df_generated.to_csv(patho_csv_path, index=False)
                    st.success("Pathogenicity data generated!")
                    st.rerun()
        else:
            patho_types = list(df_patho["Pathogenicity"].unique())
            patho_colors = {
                t: PATHO_PALETTE[i % len(PATHO_PALETTE)]
                for i, t in enumerate(patho_types)
            }

            # Graph de patho
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
                        marker_color=patho_colors.get(patho_type, "#888888"),
                    )
                )
            fig_gen_patho.update_layout(
                barmode="stack",
                barnorm="percent",
                legend=dict(yanchor="bottom", y=-1, xanchor="left", x=0.01),
                xaxis_title="Genotypes",
                yaxis_title="% Pathogenicity",
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
                    df_patho.groupby(["Cleavage Motif", "Motif Category"])
                    .size()
                    .reset_index(name="Count")
                )
                motif_counts = motif_counts[
                    motif_counts["Cleavage Motif"] != "N/A"
                ].sort_values("Count", ascending=False)
                st.dataframe(motif_counts, hide_index=True, width="stretch")

            with col_patho_2:
                st.subheader("Pathogenicity Distribution")
                vir_counts = df_patho["Pathogenicity"].value_counts()
                fig_vir_pie = go.Figure(
                    go.Pie(
                        labels=vir_counts.index,
                        values=vir_counts.values,
                        hole=0.45,
                        textposition="inside",
                        marker=dict(
                            colors=[
                                patho_colors.get(t, "#888888") for t in vir_counts.index
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

            with st.expander("Explore Pathogenecity by sequence"):
                st.dataframe(df_patho)
