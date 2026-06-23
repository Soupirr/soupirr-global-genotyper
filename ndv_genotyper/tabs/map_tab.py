"""Global Distribution Map tab."""

import os
import streamlit as st
import pandas as pd
import pydeck as pdk
import plotly.graph_objects as go
from ndv_genotyper.config import LOCATION_FOLDER, SEQ_FOLDER, PALETTE

# ============================================================================


def load_locations_database():
    df_china = pd.read_csv(os.path.join(LOCATION_FOLDER, "china_provinces.csv"))
    df_us = pd.read_csv(os.path.join(LOCATION_FOLDER, "US_capitals.csv"))
    df_world = pd.read_csv(os.path.join(LOCATION_FOLDER, "countries.csv"))
    df_russia = pd.read_csv(os.path.join(LOCATION_FOLDER, "russia_regions.csv"))
    return df_china, df_us, df_world, df_russia


df_china, df_us, df_world, df_russia = load_locations_database()

df_continent = pd.read_csv(os.path.join(LOCATION_FOLDER, "continent_map.csv"))
CONTINENT_MAP = dict(zip(df_continent["country"], df_continent["continent"]))

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


# ============================================================================


@st.cache_resource
def build_map_dataframe():
    rows = []
    for file in os.listdir(SEQ_FOLDER):
        if file.endswith(".fas"):
            path = os.path.join(SEQ_FOLDER, file)
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


# ============================================================================


def render():
    st.header("Global Distribution Map")
    st.markdown("Explore the geographic distribution of NDV genotypes worldwide.")

    if "data_loaded" not in st.session_state:
        st.info(
            "The map database contains thousands of sequences and may take a few seconds to load."
        )
        if st.button("Load Map", type="primary"):
            st.session_state["data_loaded"] = True
            st.rerun()

    if "data_loaded" in st.session_state:
        df_map = build_map_dataframe()

        # Nettoyer le DataFrame
        df_map_clean = df_map[
            (df_map["label"] != "Unknown")
            & (df_map["lat"].notna())
            & (df_map["lon"].notna())
        ].copy()

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

        df_report = df_filtered.copy()
        df_report["country"] = df_report["label"].apply(
            lambda x: x.split(",")[-1].strip()
        )
        df_report["continent"] = (
            df_report["country"].map(CONTINENT_MAP).fillna("Unknown")
        )

        df_agg = (
            df_filtered.groupby(["lat", "lon", "label"])
            .size()
            .reset_index(name="count")
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
                                [
                                    0,
                                    99,
                                    153,
                                    120,
                                ],  # rappel le 4ème nombres est l'intensité
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

        # le rapport
        if not df_filtered.empty:
            with st.expander("Genotype Distribution Report", expanded=True):
                # résumé des valeurs
                col1, col2, col3, col4 = st.columns(4)
                col1.metric("Sequences shown", len(df_report))
                col2.metric("Genotypes", df_report["genotype"].nunique())
                col3.metric("Countries", df_report["country"].nunique())
                col4.metric("Continents", df_report["continent"].nunique())

                st.divider()

                # onglets
                tab_world, tab_continent, tab_country = st.tabs(
                    ["World", "By Continent", "By Country"]
                )

                # l'onglet du monde
                with tab_world:
                    geno_counts = (
                        df_report.groupby("genotype").size().reset_index(name="count")
                    )
                    world_colors = [
                        PALETTE[i % len(PALETTE)] for i in range(len(geno_counts))
                    ]
                    world_key = "_".join(geno_counts["genotype"].tolist())

                    col_bar, col_pie = st.columns([2, 1])

                    with col_bar:
                        fig_bar = go.Figure(
                            go.Bar(
                                x=geno_counts["genotype"],
                                y=geno_counts["count"],
                                marker=dict(color=world_colors),
                                text=geno_counts["count"],
                                textposition="inside",
                            )
                        )
                        fig_bar.update_layout(
                            xaxis_title="Genotype",
                            yaxis_title="Sequences",
                            xaxis=dict(type="category"),
                            height=350,
                            margin=dict(l=0, r=0, t=20, b=0),
                            plot_bgcolor="rgba(0,0,0,0)",
                            paper_bgcolor="rgba(0,0,0,0)",
                        )
                        st.plotly_chart(
                            fig_bar, width="stretch", key=f"world_bar_{world_key}"
                        )

                    with col_pie:
                        fig_pie = go.Figure(
                            go.Pie(
                                labels=geno_counts["genotype"],
                                values=geno_counts["count"],
                                hole=0.4,
                                marker=dict(colors=world_colors),
                                textinfo="label+percent",
                                textposition="inside",
                                textfont=dict(color="white"),
                            )
                        )
                        fig_pie.update_layout(
                            height=350,
                            margin=dict(l=0, r=0, t=20, b=0),
                            paper_bgcolor="rgba(0,0,0,0)",
                            font=dict(color="white"),
                            showlegend=False,
                        )
                        st.plotly_chart(
                            fig_pie, width="stretch", key=f"world_pie_{world_key}"
                        )

                # onglet par continent
                with tab_continent:
                    cont_geno = (
                        df_report.groupby(["continent", "genotype"])
                        .size()
                        .reset_index(name="count")
                    )

                    all_genotypes = sorted(cont_geno["genotype"].unique())
                    geno_color = {
                        g: PALETTE[i % len(PALETTE)]
                        for i, g in enumerate(all_genotypes)
                    }

                    for continent in sorted(cont_geno["continent"].unique()):
                        df_c = cont_geno[cont_geno["continent"] == continent].copy()
                        df_c["pct"] = (df_c["count"] / df_c["count"].sum() * 100).round(
                            1
                        )
                        st.subheader(continent)
                        col_chart, col_stats = st.columns([3, 1])

                        with col_chart:
                            fig_cont = go.Figure(
                                go.Bar(
                                    x=df_c["genotype"],
                                    y=df_c["count"],
                                    marker=dict(
                                        color=[geno_color[g] for g in df_c["genotype"]],
                                    ),
                                    text=df_c["count"],
                                    textposition="inside",
                                )
                            )
                            fig_cont.update_layout(
                                xaxis_title="Genotype",
                                yaxis_title="Sequences",
                                xaxis=dict(type="category"),
                                height=300,
                                margin=dict(l=0, r=0, t=10, b=0),
                                plot_bgcolor="rgba(0,0,0,0)",
                                paper_bgcolor="rgba(0,0,0,0)",
                            )
                            st.plotly_chart(
                                fig_cont, width="stretch", key=f"fig_cont_{continent}"
                            )

                        with col_stats:
                            st.dataframe(
                                df_c[["genotype", "count", "pct"]]
                                .rename(
                                    columns={
                                        "genotype": "Genotype",
                                        "count": "Sequences",
                                        "pct": "%",
                                    }
                                )
                                .sort_values("%", ascending=False),
                                hide_index=True,
                                width="stretch",
                                height=300,
                            )

                        st.divider()

                # onglet par pays
                with tab_country:
                    selected_continent = st.selectbox(
                        "Select a continent",
                        sorted(df_report["continent"].unique()),
                    )
                    df_country = df_report[df_report["continent"] == selected_continent]
                    country_geno = (
                        df_country.groupby(["country", "genotype"])
                        .size()
                        .reset_index(name="count")
                    )

                    fig_country = go.Figure()
                    for geno in country_geno["genotype"].unique():
                        df_g = country_geno[country_geno["genotype"] == geno]
                        fig_country.add_trace(
                            go.Bar(
                                name=geno,
                                x=df_g["country"],
                                y=df_g["count"],
                                hovertemplate="<b>%{x}</b><br>Genotype: %{fullData.name}<br>Sequences: %{y}<extra></extra>",
                            )
                        )
                    fig_country.update_layout(
                        barmode="stack",
                        xaxis_title="Country",
                        yaxis_title="Sequences",
                        height=400,
                        margin=dict(l=0, r=0, t=20, b=0),
                        plot_bgcolor="rgba(0,0,0,0)",
                        paper_bgcolor="rgba(0,0,0,0)",
                        xaxis_tickangle=-45,
                        hoverlabel=dict(
                            bgcolor="#0c1a24",
                            font_color="white",
                            bordercolor="rgba(255,255,255,0.2)",
                        ),
                    )
                    st.plotly_chart(fig_country, width="stretch")
