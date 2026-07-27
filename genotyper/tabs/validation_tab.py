"""Validation tab : Repeats N times: randomly holds out sequences from entry's own ref dataset and tests them again's what's left"""

import os
import random
import re
import time

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from genotyper.analyzer import GenotypeIdentifier
from genotyper.tabs.analyze_tab import load_all_references

RESULT_PREFIX = "validation_results_"


def _stratified_holdout(
    references, holdout_size, min_remaining_pct=0.20, well_covered_threshold=15
):
    """Tirage aléatoire stratifié :
    - Au moins 1 séquence par génotype bien couvert (>= well_covered_threshold) est garantie
    - Au moins min_remaining_pct des séquences de chaque génotype restent dans le pool
    """
    by_genotype = {}
    for h in references:
        parts = h.split("|")
        geno = parts[2] if len(parts) >= 3 else "Unknown"
        by_genotype.setdefault(geno, []).append(h)

    guaranteed = []
    remaining_pool = {geno: list(headers) for geno, headers in by_genotype.items()}

    # Étape 1 : garantir au moins 1 séquence par génotype bien couvert
    for geno, headers in by_genotype.items():
        max_holdout = max(
            0, len(headers) - max(1, int(len(headers) * min_remaining_pct))
        )
        if len(headers) >= well_covered_threshold and max_holdout > 0:
            chosen = random.choice(remaining_pool[geno])
            guaranteed.append(chosen)
            remaining_pool[geno].remove(chosen)

    # Si on a déjà atteint ou dépassé holdout_size avec les garantis
    if len(guaranteed) >= holdout_size:
        random.shuffle(guaranteed)
        return guaranteed[:holdout_size]

    # Étape 2 : remplir le reste aléatoirement dans les slots disponibles
    extra_slots = holdout_size - len(guaranteed)
    candidates = []
    for geno, headers in remaining_pool.items():
        already_taken = len(by_genotype[geno]) - len(headers)
        max_holdout = (
            max(
                0,
                len(by_genotype[geno])
                - max(1, int(len(by_genotype[geno]) * min_remaining_pct)),
            )
            - already_taken
        )
        if max_holdout > 0:
            candidates.extend(random.sample(headers, min(max_holdout, len(headers))))

    random.shuffle(candidates)
    holdout_headers = guaranteed + candidates[:extra_slots]
    random.shuffle(holdout_headers)
    return holdout_headers


def _apply_pattern(genotype, pattern):
    if not pattern:
        return genotype
    m = re.search(pattern, genotype)
    return m.group(0) if m else genotype


def run_validation(references, holdout_size, n, method, genotype_pattern=""):
    rows = []

    total_steps = n * holdout_size
    progress = st.progress(0, text="Running validation...")
    step = 0

    for run in range(n):
        holdout_headers = _stratified_holdout(references, holdout_size)
        holdout_set = set(holdout_headers)
        pool = {h: s for h, s in references.items() if h not in holdout_set}
        identifier = GenotypeIdentifier(pool)

        for header in holdout_headers:
            step += 1
            progress.progress(
                step / total_steps, text=f"Run {run + 1}/{n} - {step}/{total_steps}"
            )
            true_genotype = _apply_pattern(header.split("|")[2], genotype_pattern)
            sequence = references[header]
            matches = identifier.identify(sequence, method=method, top_n=1)
            if not matches:
                continue
            predicted_genotype = _apply_pattern(matches[0][0], genotype_pattern)
            avg_score = matches[0][1]
            rows.append(
                {
                    "Run": run + 1,
                    "Header": header,
                    "True Genotype": true_genotype,
                    "Predicted Genotype": predicted_genotype,
                    "Match": predicted_genotype == true_genotype,
                    "Score": avg_score,
                }
            )

    progress.progress(1.0, text="Validation Complete !")
    progress.empty()
    return pd.DataFrame(rows)


def render_results(df):
    if df.empty:
        st.warning("No scorable prediction...")
        return

    run_acc = df.groupby("Run")["Match"].mean() * 100
    mean_acc = run_acc.mean()
    std_acc = run_acc.std() if len(run_acc) > 1 else 0

    st.success("Validation Complete!")
    st.subheader("Accuracy")
    c1, c2, c3 = st.columns(3)
    c1.metric("Mean accuracy", f"{mean_acc:.1f}%")
    c2.metric(
        "Standard Deviation",
        f"± {std_acc:.1f}",
    )
    c3.metric("Runs", f"{len(run_acc)}", f"{len(df)} prediction total")

    fig_runs = go.Figure(
        go.Bar(
            x=[f"Run {r}" for r in run_acc.index],
            y=run_acc.values,
            marker_color="#00c9a7",
        )
    )
    fig_runs.update_layout(
        yaxis_title="Accuracy (%)",
        height=280,
        margin=dict(l=0, r=0, t=10, b=0),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        font=dict(color="white"),
    )
    st.plotly_chart(fig_runs, width="stretch")

    st.divider()

    st.subheader("Confusion Matrix (pooled)")
    genos = sorted(set(df["True Genotype"]) | set(df["Predicted Genotype"]))
    graph_cm = pd.DataFrame(0, index=genos, columns=genos)
    for _, row in df.iterrows():
        graph_cm.loc[row["True Genotype"], row["Predicted Genotype"]] += 1

    z_values = graph_cm.values.astype(float)
    n = len(genos)
    diag_mask = np.eye(n, dtype=bool)

    z_diag = np.where(diag_mask, z_values, np.nan)
    z_off = np.where(~diag_mask, z_values, np.nan)
    text_diag = np.where(diag_mask, graph_cm.values, "").tolist()
    text_off = np.where(~diag_mask, graph_cm.values, "").tolist()

    max_off = float(np.nanmax(z_off)) if not np.all(np.isnan(z_off)) else 1.0
    max_diag = float(np.nanmax(z_diag)) if not np.all(np.isnan(z_diag)) else 1.0

    fig_cm = go.Figure()
    fig_cm.add_trace(
        go.Heatmap(
            z=z_off,
            x=genos,
            y=genos,
            colorscale=[[0, "#ffffff"], [1, "#c0392b"]],
            showscale=False,
            text=text_off,
            texttemplate="%{text}",
            textfont={"size": 11, "color": "#111111"},
            zmin=0,
            zmax=max_off,
            customdata=graph_cm.values,
            hovertemplate="True: %{y}<br>Predicted: %{x}<br>Occurrences: %{customdata}<extra></extra>",
        )
    )
    fig_cm.add_trace(
        go.Heatmap(
            z=z_diag,
            x=genos,
            y=genos,
            colorscale=[[0, "#b2ede5"], [1, "#00c9a7"]],
            showscale=False,
            text=text_diag,
            texttemplate="%{text}",
            textfont={"size": 11, "color": "#111111"},
            zmin=0,
            zmax=max_diag,
            customdata=graph_cm.values,
            hovertemplate="True: %{y}<br>Predicted: %{x}<br>Occurrences: %{customdata}<extra></extra>",
        )
    )
    fig_cm.update_layout(
        xaxis_title="Predicted",
        yaxis_title="True",
        height=max(420, len(genos) * 26),
        margin=dict(l=0, r=0, t=0, b=0),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        font=dict(color="white"),
        xaxis=dict(type="category"),
        yaxis=dict(type="category", autorange="reversed"),
    )
    st.plotly_chart(fig_cm, width="stretch")

    st.divider()

    st.subheader("Most Frequent Confustion")
    df_miss = df[~df["Match"]]
    if df_miss.empty:
        st.success("No misclassification across any run.")
    else:
        top_confusion = (
            df_miss.groupby(["True Genotype", "Predicted Genotype"])
            .size()
            .reset_index(name="Count")
            .sort_values("Count", ascending=False)
            .head(20)
        )
        st.dataframe(top_confusion, width="stretch", hide_index=True)

    st.divider()

    st.subheader("Export")
    if st.button("Save validation report", type="primary"):
        import sys

        if getattr(sys, "frozen", False):
            export_dir = os.path.join(
                os.path.dirname(sys.executable), "exports", "validation"
            )
        else:
            export_dir = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "..",
                "..",
                "exports",
                "validation",
            )
        os.makedirs(export_dir, exist_ok=True)

        timestamp = time.strftime("%Y%m%d_%H%M%S")

        # Fichier 1 : prédictions brutes
        df.to_csv(
            os.path.join(export_dir, f"validation_predictions_{timestamp}.csv"),
            index=False,
        )

        # Fichier 2 : résumé global
        run_acc = df.groupby("Run")["Match"].mean() * 100
        mean_acc = run_acc.mean()
        std_acc = run_acc.std() if len(run_acc) > 1 else 0

        df_summary = pd.DataFrame(
            {
                "Metric": [
                    "Mean Accuracy (%)",
                    "Standard Deviation",
                    "Number of Runs",
                    "Total Predictions",
                ],
                "Value": [f"{mean_acc:.2f}", f"{std_acc:.2f}", len(run_acc), len(df)],
            }
        )
        df_per_run = run_acc.reset_index()
        df_per_run.columns = ["Run", "Accuracy (%)"]
        df_per_run["Accuracy (%)"] = df_per_run["Accuracy (%)"].map("{:.2f}".format)

        df_miss = df[~df["Match"]]
        top_confusion = (
            (
                df_miss.groupby(["True Genotype", "Predicted Genotype"])
                .size()
                .reset_index(name="Count")
                .sort_values("Count", ascending=False)
                .head(20)
            )
            if not df_miss.empty
            else pd.DataFrame(columns=["True Genotype", "Predicted Genotype", "Count"])
        )

        summary_path = os.path.join(export_dir, f"validation_summary_{timestamp}.csv")
        with open(summary_path, "w", encoding="utf-8") as f:
            f.write("=== GLOBAL METRICS ===\n")
            df_summary.to_csv(f, index=False)
            f.write("\n=== ACCURACY PER RUN ===\n")
            df_per_run.to_csv(f, index=False)
            f.write("\n=== TOP CONFUSIONS ===\n")
            top_confusion.to_csv(f, index=False)

        st.success(f"Saved to: {export_dir}")


def render(path, entry_config=None):
    st.header("Precision Validation")
    st.markdown(
        "Randomly holds out sequences from this entry's own reference dataset "
        "and tests the genotyper against them, repeated multiple times to get "
        "a mean accuracy instead of a single lucky/unlucky draw."
    )

    db_references, _, total_count, _ = load_all_references(path)
    is_multi = bool(
        db_references and isinstance(next(iter(db_references.values())), dict)
    )

    if is_multi:
        gene_names = list(db_references.keys())
        selected_gene = st.segmented_control("Gene", gene_names, default=gene_names[0])
        if selected_gene not in gene_names:
            selected_gene = gene_names[0]
        references = db_references[selected_gene]
        gene_total = len(references)
        result_key = f"{RESULT_PREFIX}{os.path.basename(path)}_{selected_gene}"
        gene_cfg = (
            (entry_config.get("genes", {}) or {}).get(selected_gene, {})
            if entry_config
            else {}
        )
        genotype_pattern = gene_cfg.get("genotype_pattern", "")
    else:
        references = db_references
        gene_total = total_count
        result_key = f"{RESULT_PREFIX}{os.path.basename(path)}"
        genotype_pattern = ""

    if gene_total < 10:
        st.warning(
            f"Only {gene_total} reference sequences in this entry, too few for a meaningful holdout test."
        )
        return

    st.info(f"Reference dataset : **{gene_total}** sequences")

    col1, col2 = st.columns(2)
    with col1:
        default_holdout = max(1, min(100, gene_total // 3))
        holdout_size = st.slider(
            "Sequences held out per run",
            min_value=1,
            max_value=max(1, gene_total // 2),
            value=default_holdout,
            help="Removed from the matching pool each run, then tested against what's left.",
        )
        n = st.number_input(
            "Number of runs", min_value=1, max_value=100, value=5, step=1
        )

    with col2:
        method_label = st.radio("Similarity Method", ["Pairwise", "Hamming"])
        if "Pairwise" in method_label:
            st.caption("Pairwise is really slow ~1min per sequences.")
        if "Hamming" in method_label:
            st.caption(
                "Only use hamming if you are sure that both the reference sequences and the sequences you will analyse are aligned between themselves."
            )

    run = st.button("Run Validation Test", type="primary")
    clear = st.button("Clear results")
    if clear:
        st.session_state.pop(result_key, None)

    if run:
        method = "hamming" if "Hamming" in method_label else "pairwise"
        df = run_validation(references, holdout_size, int(n), method, genotype_pattern)
        st.session_state[result_key] = df
    if result_key in st.session_state:
        st.divider()
        render_results(st.session_state[result_key])
