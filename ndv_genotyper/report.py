"""
Assembles the analysis results, the similarity matrix, and the phylogenetic
tree(s) into one downloadable HTML file.
"""

import html as _html
from datetime import datetime
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs


# Dark theme matching the app (so embedded charts blend in).
_CSS = """
  body { background:#060d14; color:#b8d4e8; font-family:'Segoe UI',system-ui,sans-serif;
         margin:0; padding:32px 40px; }
  h1 { color:#00c9a7; font-weight:300; letter-spacing:2px; }
  h2 { color:#0099cc; font-weight:300; margin-top:36px;
       border-bottom:1px solid #0e2535; padding-bottom:6px; }
  h3 { color:#6699cc; font-weight:400; }
  .meta { color:#5e7a8f; font-size:0.9em; }
  .card { background:#0c1a24; border:1px solid #0e2535; border-radius:10px;
          padding:18px 22px; margin:16px 0; }
  table.rpt { border-collapse:collapse; width:100%; margin:10px 0; font-size:0.9em; }
  table.rpt th { background:#0c1a24; color:#00c9a7; text-align:left;
                 padding:8px 12px; border-bottom:1px solid #0e2535; font-weight:500; }
  table.rpt td { padding:7px 12px; border-bottom:1px solid #0e2535; }
  .patho-virulent { color:#ff6b6b; }
  .patho-low { color:#2ECC71; }
  .patho-undet { color:#888; }
"""


def _fig_matrix(fig):
    fig2 = go.Figure(fig)
    fig2.update_layout(
        yaxis=dict(automargin=True),
        margin=dict(l=150, r=50, t=40, b=100),
    )
    return fig2.to_html(
        full_html=False,
        include_plotlyjs=False,
        config={"displayModeBar": True, "displaylogo": False, "scrollZoom": False},
    )


def _fig_tree(fig):
    fig2 = go.Figure(fig)
    fig2.update_layout(width=None, autosize=True)
    return fig2.to_html(
        full_html=False,
        include_plotlyjs=False,
        config={"displayModeBar": True, "displaylogo": False, "scrollZoom": True},
    )


def _table(df):
    return df.to_html(index=False, border=0, classes="rpt", escape=True)


def build_report_html(
    all_results,
    export_df,
    method,
    elapsed_time,
    matrix_fig=None,
    tree_figs=None,
):
    tree_figs = tree_figs or {}
    has_figs = matrix_fig is not None or bool(tree_figs)

    p = [
        "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        "<title>NDV Genotyper Report</title>",
    ]
    if has_figs:
        p.append(f"<script>{get_plotlyjs()}</script>")
    p.append(f"<style>{_CSS}</style></head><body>")

    # --- header ---
    p.append("<h1>NDV Genotyper — Analysis Report</h1>")
    p.append(
        f"<p class='meta'>Generated {datetime.now():%Y-%m-%d %H:%M} &middot; "
        f"Method: {_html.escape(str(method)).title()} &middot; "
        f"{len(all_results)} sequence(s) &middot; {elapsed_time:.2f}s</p>"
    )

    # --- full analysis details table ---
    p.append("<h2>Full Analysis Details</h2>")
    p.append(_table(export_df))

    # --- similarity matrix ---
    if matrix_fig is not None:
        p.append("<h2>Similarity Matrix</h2>")
        p.append(_fig_matrix(matrix_fig))

    # --- phylogenetic tree(s) ---
    if tree_figs:
        p.append("<h2>Phylogenetic Tree(s)</h2>")
        for title, fig in tree_figs.items():
            p.append(f"<h3>{_html.escape(str(title))}</h3>")
            p.append(_fig_tree(fig))
    else:
        p.append("<h2>Phylogenetic Tree(s)</h2>")
        p.append(
            "<p class='meta'>No trees were generated. Build them in the "
            "Phylogenetic Tree tab before generating the report to include them.</p>"
        )

    p.append(
        "<p class='meta' style='margin-top:40px;'>Pathogenicity predictions are "
        "motif-based (Wang et al., 2017) and indicative only — validate by assay "
        "before any clinical or regulatory conclusion.</p>"
    )
    p.append("</body></html>")
    return "".join(p)
