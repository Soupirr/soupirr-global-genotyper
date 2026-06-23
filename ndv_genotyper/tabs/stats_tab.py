"""Statistics / Help tab."""

import os
import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from ndv_genotyper.tabs import map_tab, analyze_tab
from ndv_genotyper.config import DATA_FOLDER, PATHO_FOLDER, MISC_FOLDER, PALETTE

# ============================================================================


def render():
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
        with open(
            os.path.join(MISC_FOLDER, "QUICK_START.md"), "r", encoding="utf-8"
        ) as f:
            informations = f.read()

        st.markdown(informations)

    with path_tab:
        with open(
            os.path.join(MISC_FOLDER, "PATHOGENICITY_CRITERIA.md"),
            "r",
            encoding="utf-8",
        ) as p:
            pathogénicité = p.read()

        st.markdown(pathogénicité)

    with stat_tab:
        st.header("Database Statistics")

        if "data_loaded" not in st.session_state:
            st.info(
                "Building the statistics scans every sequence in the database. "
                "Click to load."
            )
            if st.button("Load Statistics", type="primary"):
                st.session_state["data_loaded"] = True
                st.rerun()
            return

        db_references, db_files_count, db_total_count, _ = (
            analyze_tab.load_all_references()
        )

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
