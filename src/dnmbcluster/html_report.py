"""Auto-generated HTML summary report.

Consolidates all pipeline outputs — key numbers, plot thumbnails,
and download links — into a single browsable HTML file. No external
template engine needed; pure stdlib string formatting.
"""
from __future__ import annotations

import base64
import logging
from pathlib import Path

import pyarrow.parquet as pq

log = logging.getLogger(__name__)


def _img_tag(pdf_path: Path, width: int = 600) -> str:
    """Embed a PDF as a base64 <object> or link if the file exists."""
    if not pdf_path.exists():
        return ""
    return (
        f'<div class="plot-card">'
        f'<h3>{pdf_path.stem}</h3>'
        f'<object data="{pdf_path.name}" type="application/pdf" '
        f'width="{width}" height="{int(width * 0.75)}">'
        f'<a href="{pdf_path.name}">{pdf_path.stem} (PDF)</a>'
        f'</object></div>'
    )


def _stat_row(label: str, value) -> str:
    return f"<tr><td>{label}</td><td><strong>{value}</strong></td></tr>"


def write_html_report(results_dir: Path) -> Path:
    """Generate ``results_dir/report.html``."""
    processed = results_dir / "dnmb" / "processed"
    inputs = results_dir / "dnmb" / "inputs"
    plots_dir = results_dir / "plots"
    out_path = results_dir / "report.html"

    # --- Collect stats ---
    stats_rows = []
    meta = None
    if (inputs / "genome_meta.parquet").exists():
        meta = pq.read_table(inputs / "genome_meta.parquet")
        stats_rows.append(_stat_row("Genomes", f"{meta.num_rows:,}"))

    if (processed / "clusters.parquet").exists():
        cl = pq.read_table(processed / "clusters.parquet")
        stats_rows.append(_stat_row("Total CDS", f"{cl.num_rows:,}"))

    if (processed / "cluster_summary.parquet").exists():
        sm = pq.read_table(processed / "cluster_summary.parquet").to_pydict()
        cats = sm.get("category", [])
        n_cl = len(cats)
        n_core = cats.count("core")
        n_acc = cats.count("accessory")
        n_uniq = cats.count("unique")
        stats_rows.append(_stat_row("Clusters", f"{n_cl:,}"))
        stats_rows.append(_stat_row("Core", f"{n_core:,} ({n_core/max(n_cl,1)*100:.1f}%)"))
        stats_rows.append(_stat_row("Accessory", f"{n_acc:,} ({n_acc/max(n_cl,1)*100:.1f}%)"))
        stats_rows.append(_stat_row("Unique", f"{n_uniq:,} ({n_uniq/max(n_cl,1)*100:.1f}%)"))

    if (processed / "phylo_tree.nwk").exists():
        stats_rows.append(_stat_row("Phylo tree", "yes"))

    if (processed / "ani_matrix.parquet").exists():
        ani = pq.read_table(processed / "ani_matrix.parquet")
        vals = [v for v in ani.column("ani_percent").to_pylist() if v < 99.99]
        if vals:
            stats_rows.append(_stat_row("ANI range", f"{min(vals):.1f}% – {max(vals):.1f}%"))

    # --- Plot thumbnails ---
    plot_cards = []
    if plots_dir.exists():
        for pdf in sorted(plots_dir.glob("*.pdf")):
            plot_cards.append(
                f'<div class="plot-card">'
                f'<h3>{pdf.stem.replace("_", " ").title()}</h3>'
                f'<a href="plots/{pdf.name}" target="_blank">'
                f'<div class="plot-placeholder">📊 {pdf.name}</div>'
                f'</a></div>'
            )

    # --- Build HTML ---
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>DNMBcluster Report</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         max-width: 1200px; margin: 0 auto; padding: 20px; background: #fafafa; }}
  h1 {{ color: #2C5F7A; border-bottom: 2px solid #2C5F7A; padding-bottom: 8px; }}
  h2 {{ color: #3C5488; margin-top: 30px; }}
  table {{ border-collapse: collapse; margin: 10px 0; }}
  td {{ padding: 6px 16px; border-bottom: 1px solid #eee; }}
  td:first-child {{ color: #606060; }}
  .plots-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(350px, 1fr));
                 gap: 16px; margin-top: 16px; }}
  .plot-card {{ background: white; border-radius: 8px; padding: 12px;
               box-shadow: 0 1px 4px rgba(0,0,0,0.1); }}
  .plot-card h3 {{ margin: 0 0 8px; font-size: 14px; color: #3C5488; }}
  .plot-placeholder {{ background: #f0f0f0; border-radius: 4px; padding: 40px 20px;
                       text-align: center; color: #808080; font-size: 13px; }}
  .plot-placeholder:hover {{ background: #e8e8e8; }}
  a {{ color: #2C5F7A; text-decoration: none; }}
  a:hover {{ text-decoration: underline; }}
  .footer {{ margin-top: 40px; padding-top: 12px; border-top: 1px solid #ddd;
             color: #999; font-size: 12px; }}
</style>
</head>
<body>
<h1>DNMBcluster Report</h1>
<p><strong>Results directory:</strong> <code>{results_dir}</code></p>

<h2>Summary Statistics</h2>
<table>{"".join(stats_rows)}</table>

<h2>Visualizations ({len(plot_cards)} plots)</h2>
<div class="plots-grid">{"".join(plot_cards)}</div>

<h2>Output Files</h2>
<ul>
  <li><a href="gene_presence_absence.csv">gene_presence_absence.csv</a> (Roary-compatible)</li>
  <li><a href="gene_presence_absence.Rtab">gene_presence_absence.Rtab</a></li>
  <li><a href="comparative_genomics_10.xlsx">comparative_genomics_*.xlsx</a> (merged sheet)</li>
</ul>

<div class="footer">
  Generated by DNMBcluster · <a href="https://github.com/JAEYOONSUNG/DNMBcluster">GitHub</a>
</div>
</body>
</html>"""

    out_path.write_text(html, encoding="utf-8")
    log.info("HTML report written to %s", out_path)
    return out_path
