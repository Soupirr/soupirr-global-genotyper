"""Validation / accuracy testing tab (hidden dev tool).

Runs the genotyper over a labeled holdout set and reports:
  - hierarchical accuracy (exact sub-genotype / top genotype / class)
  - micro vs macro accuracy (honest under class imbalance)
  - per-class (I/II) breakdown
  - a misclassifications-only inspector
  - confusion matrix + per-genotype precision/recall/F1
  - top-match similarity distribution for correct vs incorrect calls

Results are cached in session_state so changing a view doesn't re-run the analysis.
"""

import os
import random
from datetime import datetime

import streamlit as st
import pandas as pd
import plotly.graph_objects as go

from ndv_genotyper.analyzer import FASTAParser, analyze_newcastle_sequence
from ndv_genotyper.config import VALIDATION_FOLDER

TEST_FASTA = os.path.join(VALIDATION_FOLDER, "sequences_cleaned.fasta")
TEST_CSV = os.path.join(VALIDATION_FOLDER, "genotypes_cleaned.csv")

RESULT_KEY = "validation_results"


# ---------------------------------------------------------------------------
# Genotype hierarchy helpers
# Class I = Arabic leading digit (1.x), Class II = Roman (matches analyzer.get_class)
# ---------------------------------------------------------------------------
def genotype_class(g):
    return "Class I" if str(g)[:1].isdigit() else "Class II"


def top_genotype(g):
    """Top-level genotype: 'VII.1.1' -> 'VII', '1.1.2' -> '1'."""
    return str(g).split(".")[0]


def match_level(true, pred):
    """How close a prediction is, on the nested NDV nomenclature."""
    if pred == true:
        return "Exact"
    if top_genotype(pred) == top_genotype(true):
        return "Genotype"  # same top-level genotype, wrong sub-lineage
    if genotype_class(pred) == genotype_class(true):
        return "Class"  # only the class is right
    return "Wrong"


@st.cache_data
def load_test_data():
    """Load test sequences and expected genotypes."""
    try:
        df_test = pd.read_csv(TEST_CSV)
        test_genotypes = dict(zip(df_test["ID"], df_test["Genotype"]))
        test_seqs = FASTAParser.parse_file(TEST_FASTA)
        return test_seqs, test_genotypes, len(df_test)
    except FileNotFoundError:
        return {}, {}, 0


# ---------------------------------------------------------------------------
# Core run
# ---------------------------------------------------------------------------
def run_validation(
    test_sequences, test_genotypes, references, selected_ids, method, loo
):
    progress_bar = st.progress(0, text="Starting validation...")
    predictions = []

    for idx, seq_id in enumerate(selected_ids):
        progress_bar.progress(
            (idx + 1) / len(selected_ids),
            text=f"Analyzing {idx + 1}/{len(selected_ids)}: {seq_id}",
        )

        sequence = test_sequences[seq_id]
        true_genotype = test_genotypes.get(seq_id, "Unknown")
        if true_genotype == "Unknown":
            continue  # no label -> can't score, skip

        # leave-one-out: drop any reference identical to the query
        refs = (
            {h: s for h, s in references.items() if s != sequence}
            if loo
            else references
        )

        single_fasta = f">{seq_id}\n{sequence}"
        try:
            res = analyze_newcastle_sequence(
                input_fasta=single_fasta,
                reference_sequences=refs,
                top_matches=1,
                similarity_method=method,
            )
            if res["genotype_matches"]:
                pred = res["genotype_matches"][0][0]
                predictions.append(
                    {
                        "ID": seq_id,
                        "True Genotype": true_genotype,
                        "Predicted Genotype": pred,
                        "Match": pred == true_genotype,
                        "Level": match_level(true_genotype, pred),
                        "True Class": genotype_class(true_genotype),
                        "Similarity": res["genotype_matches"][0][1],
                    }
                )
        except Exception as e:
            st.error(f"Error analyzing {seq_id}: {e}")

    progress_bar.progress(1.0, text="Analysis complete!")
    progress_bar.empty()
    return pd.DataFrame(predictions)


# ---------------------------------------------------------------------------
# Results rendering (reads from session_state, so it survives unrelated reruns)
# ---------------------------------------------------------------------------
def render_results(df):
    total = len(df)
    if total == 0:
        st.warning("No scorable predictions were produced.")
        return

    levels = df["Level"]
    exact = (levels == "Exact").sum()
    geno_ok = levels.isin(["Exact", "Genotype"]).sum()
    class_ok = levels.isin(["Exact", "Genotype", "Class"]).sum()

    # per-genotype recall (for macro average)
    recalls = []
    for g in df["True Genotype"].unique():
        sub = df[df["True Genotype"] == g]
        recalls.append((sub["Match"].sum() / len(sub)) * 100)
    macro = sum(recalls) / len(recalls) if recalls else 0

    st.success("✅ Validation complete!")

    # --- A + B: hierarchical + macro accuracy ---
    st.subheader("Accuracy")
    c1, c2, c3 = st.columns(3)
    c1.metric(
        "Exact (sub-genotype)",
        f"{exact / total * 100:.1f}%",
        f"{exact}/{total}",
        help="Predicted sub-genotype matches exactly (e.g. VII.1.1 == VII.1.1)",
    )
    c2.metric(
        "Top genotype",
        f"{geno_ok / total * 100:.1f}%",
        help="Correct at the top-level genotype (e.g. VII.* counted as VII)",
    )
    c3.metric(
        "Macro accuracy",
        f"{macro:.1f}%",
        help="Mean of per-genotype recall — unweighted, so rare genotypes count "
        "as much as common ones. Compare with Exact (micro) to see imbalance impact.",
    )

    # --- C: per-class breakdown ---
    st.subheader("By Class")
    cols = st.columns(2)
    for col, cls in zip(cols, ["Class I", "Class II"]):
        sub = df[df["True Class"] == cls]
        if len(sub):
            acc = sub["Match"].mean() * 100
            col.metric(
                f"{cls} (exact)", f"{acc:.1f}%", f"{int(sub['Match'].sum())}/{len(sub)}"
            )
        else:
            col.metric(f"{cls} (exact)", "—", "no samples")

    st.divider()

    # --- F: misclassifications-only inspector ---
    st.subheader("Misclassifications")
    df_miss = df[~df["Match"]].copy()
    if df_miss.empty:
        st.success("No misclassifications 🎉")
    else:
        df_miss = df_miss.rename(columns={"Level": "Closeness"})
        view = df_miss[
            ["ID", "True Genotype", "Predicted Genotype", "Closeness", "Similarity"]
        ].sort_values("Closeness")

        level_colors = {
            "Genotype": "background-color:#8a6d0b",  # right genotype, wrong sub
            "Class": "background-color:#9c4221",  # only class right
            "Wrong": "background-color:#7a1730",  # fully wrong
        }
        styled = view.style.map(
            lambda v: level_colors.get(v, ""), subset=["Closeness"]
        ).format({"Similarity": "{:.2f}"})
        st.dataframe(styled, width="stretch", hide_index=True)
        st.caption(
            "Closeness = how near the miss was: **Genotype** (right top-level, wrong "
            "sub-lineage) · **Class** (only I/II correct) · **Wrong** (different class)."
        )

    st.divider()

    # --- E: similarity distribution, correct vs incorrect ---
    st.subheader("Top-match Similarity: correct vs incorrect")
    fig_sim = go.Figure()
    fig_sim.add_trace(
        go.Histogram(
            x=df[df["Match"]]["Similarity"],
            name="Correct",
            marker_color="#2ECC71",
            opacity=0.7,
            xbins=dict(size=1),
        )
    )
    fig_sim.add_trace(
        go.Histogram(
            x=df[~df["Match"]]["Similarity"],
            name="Incorrect",
            marker_color="#FF4444",
            opacity=0.7,
            xbins=dict(size=1),
        )
    )
    fig_sim.update_layout(
        barmode="overlay",
        xaxis_title="Top-match similarity (%)",
        yaxis_title="Count",
        height=320,
        margin=dict(l=0, r=0, t=10, b=0),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        font=dict(color="white"),
        legend=dict(orientation="h", y=1.1),
    )
    st.plotly_chart(fig_sim, width="stretch")
    st.caption(
        "If incorrect calls cluster at lower similarity, a confidence threshold "
        "(e.g. the app's 95% warning) is justified."
    )

    st.divider()

    # --- confusion matrix (sorted by class, then genotype) ---
    st.subheader("Confusion Matrix")
    genos = sorted(
        set(df["True Genotype"]) | set(df["Predicted Genotype"]),
        key=lambda g: (genotype_class(g), top_genotype(g), str(g)),
    )
    cm = pd.DataFrame(0, index=genos, columns=genos)
    for _, row in df.iterrows():
        cm.loc[row["True Genotype"], row["Predicted Genotype"]] += 1

    fig_cm = go.Figure(
        data=go.Heatmap(
            z=cm.values,
            x=cm.columns,
            y=cm.index,
            colorscale="Blues",
            text=cm.values,
            texttemplate="%{text}",
            textfont={"size": 11},
        )
    )
    fig_cm.update_layout(
        xaxis_title="Predicted",
        yaxis_title="True",
        height=max(420, len(genos) * 26),
        margin=dict(l=0, r=0, t=10, b=0),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        font=dict(color="white"),
        yaxis=dict(autorange="reversed"),
    )
    st.plotly_chart(fig_cm, width="stretch")

    st.divider()

    # --- per-genotype precision / recall / F1 ---
    st.subheader("Per-Genotype Metrics")
    rows = []
    for g in genos:
        mask_true = df["True Genotype"] == g
        mask_pred = df["Predicted Genotype"] == g
        mask_correct = mask_true & df["Match"]
        n_true, n_pred = mask_true.sum(), mask_pred.sum()
        recall = (mask_correct.sum() / n_true * 100) if n_true else 0
        precision = (mask_correct.sum() / n_pred * 100) if n_pred else 0
        f1 = (
            2 * precision * recall / (precision + recall) if (precision + recall) else 0
        )
        if n_true or n_pred:
            rows.append(
                {
                    "Genotype": g,
                    "Class": genotype_class(g),
                    "Samples": int(n_true),
                    "Precision": f"{precision:.1f}%",
                    "Recall": f"{recall:.1f}%",
                    "F1-Score": f"{f1:.1f}%",
                }
            )
    st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)

    st.divider()

    # --- G: export (predictions + summary header) ---
    st.subheader("Export")
    summary = pd.DataFrame(
        [
            {"Metric": "Exact accuracy", "Value": f"{exact / total * 100:.1f}%"},
            {
                "Metric": "Top-genotype accuracy",
                "Value": f"{geno_ok / total * 100:.1f}%",
            },
            {"Metric": "Class accuracy", "Value": f"{class_ok / total * 100:.1f}%"},
            {"Metric": "Macro accuracy", "Value": f"{macro:.1f}%"},
            {"Metric": "N tested", "Value": str(total)},
        ]
    )
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    e1, e2 = st.columns(2)
    e1.download_button(
        "Download predictions (CSV)",
        data=df.to_csv(index=False),
        file_name=f"validation_predictions_{ts}.csv",
        mime="text/csv",
        width="stretch",
    )
    e2.download_button(
        "Download summary (CSV)",
        data=summary.to_csv(index=False),
        file_name=f"validation_summary_{ts}.csv",
        mime="text/csv",
        width="stretch",
    )


# ---------------------------------------------------------------------------
# Tab entry point
# ---------------------------------------------------------------------------
def render(references):
    st.header("Test Performance & Validation")
    st.markdown(
        "Validate the genotyper against a labeled holdout set and measure performance."
    )

    test_sequences, test_genotypes, total_test = load_test_data()

    if not test_sequences:
        st.warning(
            "⚠️ Test files not found. Expected `sequences_cleaned.fasta` and "
            "`genotypes_cleaned.csv` in `data/validation/`."
        )
        return

    # dataset overview (with class split)
    classes = pd.Series([genotype_class(g) for g in test_genotypes.values()])
    n_c1 = int((classes == "Class I").sum())
    n_c2 = int((classes == "Class II").sum())
    st.info(
        f"Holdout dataset: **{len(test_sequences)}** sequences "
        f"({n_c2} Class II, {n_c1} Class I) · "
        f"{len(set(test_genotypes.values()))} genotypes"
    )

    # --- configuration ---
    col1, col2 = st.columns(2)
    with col1:
        method_label = st.radio(
            "Similarity Method",
            ["Hamming (fast)", "Pairwise (accurate)"],
            key="test_similarity",
            help="Same methods as the main analysis tab.",
        )
        loo = st.checkbox(
            "Leave-one-out",
            value=False,
            help="Exclude any reference sequence identical to the query before "
            "matching. Safety net against residual duplicates; the holdout is "
            "already leakage-free so this is optional.",
        )
    with col2:
        test_all = st.checkbox("Test entire holdout", value=True)
        sample_size = st.slider(
            "Number of sequences to test",
            min_value=1,
            max_value=len(test_sequences),
            value=len(test_sequences),
            disabled=test_all,
            help="Used only when 'Test entire holdout' is off.",
        )
        seed = st.number_input(
            "Random seed",
            value=42,
            step=1,
            disabled=test_all,
            help="Makes the sampled subset reproducible.",
        )

    run = st.button("Run Validation Test", type="primary")
    clear = st.button("Clear results")
    if clear:
        st.session_state.pop(RESULT_KEY, None)

    if run:
        method = "hamming" if "Hamming" in method_label else "pairwise"
        all_ids = list(test_sequences.keys())
        if test_all:
            selected_ids = sorted(all_ids)
        else:
            random.seed(int(seed))
            selected_ids = random.sample(all_ids, min(sample_size, len(all_ids)))

        st.info(
            f"Testing {len(selected_ids)} sequences ({method}, leave-one-out={loo})..."
        )
        df = run_validation(
            test_sequences, test_genotypes, references, selected_ids, method, loo
        )
        st.session_state[RESULT_KEY] = df

    if RESULT_KEY in st.session_state:
        st.divider()
        render_results(st.session_state[RESULT_KEY])
