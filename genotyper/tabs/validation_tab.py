"""Validation tab : Repeats N times: randomly holds out sequences from entry's own ref dataset and tests them again's what's left"""

import os
import random
import streamlit as st
import pandas as pd
import plotly.graph_objects as go

from genotyper.analyzer import GenotypeIdentifier
from genotyper.tabs.analyze_tab import load_all_references

RESULT_PREFIX = "validation_results_"


def run_validation(references, holdout_size, n, method):
    all_headers = list(references.keys())
    rows = []

    total_steps = n * holdout_size
    progress = st.progress(0, text="Running validation...")
    step = 0

    for run in range(n):
        holdout_headers = random.sample(
            all_headers, min(holdout_size, len(all_headers))
        )
        holdout_set = set(holdout_headers)
        pool = {h: s for h, s in references.items() if h not in holdout_set}
        identifier = GenotypeIdentifier(pool)

        for header in holdout_headers:
            step += 1
            progress.progress(
                step / total_steps, text=f"Run {run + 1}/{n} - {step}/{total_steps}"
            )
            true_genotype = header.split("|")[2]
            sequence = references[header]
            matches = identifier.identify(sequence, method=method, top_n=1)
            if not matches:
                continue
            predicted_genotype, avg_score = matches[0][0], matches[0][1]
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

    fig_cm = go.Figure(
        data=go.Heatmap(
            z=graph_cm.values,
            x=graph_cm.columns,
            y=graph_cm.index,
            colorscale="Blues",
            text=graph_cm.values,
            texttemplate="%{text}",
            textfont={"size": 11},
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
    st.download_button(
        "Download raw prediction CSV",
        data=df.to_csv(index=False),
        file_name="validation_predictions.csv",
        mime="text/csv",
    )


def render(path):
    st.header("Precision Validation")
    st.markdown(
        "Randomly holds out sequences from this entry's own reference dataset "
        "and tests the genotyper against them, repeated multiple times to get "
        "a mean accuracy instead of a single lucky/unlucky draw."
    )

    references, _, total_count, _ = load_all_references(path)
    result_key = f"{RESULT_PREFIX}{os.path.basename(path)}"

    if total_count < 10:
        st.warning(
            f"Only {total_count} reference sequences in this entry, too few for a meaningful holdout test."
        )
        return

    st.info(f"Reference dataset : **{total_count}** sequences")

    col1, col2 = st.columns(2)
    with col1:
        default_holdout = max(1, min(100, total_count // 3))
        holdout_size = st.slider(
            "Sequences held out per run",
            min_value=1,
            max_value=max(1, total_count - 5),
            value=default_holdout,
            help="Removed from the matching pool each run, then tested against what's left.",
        )
        n = st.number_input(
            "Number of runs", min_value=1, max_value=50, value=5, step=1
        )

    with col2:
        method_label = st.radio(
            "Similarity Method", ["Pairwise (accurate)", "Hamming (fast)"]
        )
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
        df = run_validation(references, holdout_size, int(n), method)
        st.session_state[result_key] = df
    if result_key in st.session_state:
        st.divider()
        render_results(st.session_state[result_key])
