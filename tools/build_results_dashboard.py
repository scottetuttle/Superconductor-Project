"""Build a self-contained HTML browser for SHS benchmark results.

Run from the repository root:
    python tools/build_results_dashboard.py
"""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
from urllib.parse import quote


MEDIA_SUFFIXES = {".png", ".gif", ".jpg", ".jpeg", ".svg"}
DATA_SUFFIXES = {".json", ".csv"}


def _display(value):
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.6g}"
    if value is None:
        return "—"
    return str(value)


def _relative_url(path, dashboard):
    relative = path.resolve().relative_to(dashboard.parent.resolve())
    return "/".join(quote(part) for part in relative.parts)


def _metric_table(mapping):
    if not mapping:
        return '<p class="muted">No numeric summary was recorded.</p>'
    rows = []
    for key, value in mapping.items():
        if isinstance(value, (dict, list)):
            continue
        rows.append(
            f"<tr><th>{html.escape(str(key).replace('_', ' '))}</th>"
            f"<td>{html.escape(_display(value))}</td></tr>"
        )
    return '<table class="metrics">' + "".join(rows) + "</table>"


def _case_block(case, case_dir, dashboard):
    label = str(case.get("label", case_dir.name if case_dir else "case"))
    final = case.get("final", {})
    regime = case.get("model_regime", {})
    warnings = regime.get("warnings", [])
    status = "pass" if case.get("converged", False) and regime.get("passes", True) else "warn"
    files = sorted(case_dir.rglob("*")) if case_dir and case_dir.exists() else []
    media = [path for path in files if path.is_file() and path.suffix.lower() in MEDIA_SUFFIXES]
    data = [path for path in files if path.is_file() and path.suffix.lower() in DATA_SUFFIXES]
    warning_html = "".join(f"<li>{html.escape(str(item))}</li>" for item in warnings)
    gallery = "".join(
        f'<figure><a href="{_relative_url(path, dashboard)}"><img loading="lazy" '
        f'src="{_relative_url(path, dashboard)}" alt="{html.escape(path.name)}"></a>'
        f'<figcaption>{html.escape(path.name)}</figcaption></figure>'
        for path in media
    )
    links = " ".join(
        f'<a class="file" href="{_relative_url(path, dashboard)}">{html.escape(path.name)}</a>'
        for path in data
    )
    return f"""
    <details class="case" open>
      <summary><span class="badge {status}">{status.upper()}</span> {html.escape(label)}</summary>
      <div class="case-grid">
        <section><h4>Final physical and numerical values</h4>{_metric_table(final)}</section>
        <section><h4>Model regime</h4>{_metric_table(regime.get('metrics', {}))}
          {'<ul class="warnings">' + warning_html + '</ul>' if warnings else '<p class="pass-text">Configured validity checks passed.</p>'}
          <h4>Data files</h4><div class="files">{links or 'No case data files found.'}</div>
        </section>
      </div>
      <div class="gallery">{gallery or '<p class="muted">This run did not generate visual media.</p>'}</div>
    </details>"""


def build_dashboard(results_root: Path, output: Path):
    summaries = sorted(results_root.glob("*/summary.json"), key=lambda p: p.stat().st_mtime, reverse=True)
    run_blocks = []
    for summary_path in summaries:
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            run_blocks.append(f'<section class="run"><h2>{html.escape(summary_path.parent.name)}</h2><p>Unreadable summary: {html.escape(str(error))}</p></section>')
            continue
        cases = summary.get("cases", [])
        case_blocks = []
        for case in cases:
            case_dir = summary_path.parent / str(case.get("label", ""))
            case_blocks.append(_case_block(case, case_dir, output))
        root_media = [path for path in summary_path.parent.iterdir() if path.is_file() and path.suffix.lower() in MEDIA_SUFFIXES]
        root_gallery = "".join(
            f'<figure><a href="{_relative_url(path, output)}"><img loading="lazy" src="{_relative_url(path, output)}" alt="{html.escape(path.name)}"></a><figcaption>{html.escape(path.name)}</figcaption></figure>'
            for path in root_media
        )
        run_blocks.append(f"""
        <section class="run">
          <h2>{html.escape(summary_path.parent.name)}</h2>
          <p><a class="file" href="{_relative_url(summary_path, output)}">combined summary.json</a> · {len(cases)} case(s)</p>
          <div class="gallery">{root_gallery}</div>
          {''.join(case_blocks) if case_blocks else '<p class="muted">No case summaries found.</p>'}
        </section>""")

    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>SHS Results Dashboard</title>
<style>
:root {{ color-scheme: dark; --bg:#0b1118; --panel:#131d28; --line:#2a3b4d; --text:#e8f0f7; --muted:#9db0c0; --accent:#63c7ff; --good:#5ee39a; --warn:#ffbf69; }}
* {{ box-sizing:border-box }} body {{ margin:0; background:var(--bg); color:var(--text); font:15px/1.45 system-ui,sans-serif }}
header,main {{ max-width:1500px; margin:auto; padding:24px }} header h1 {{ margin-bottom:4px }} .muted {{ color:var(--muted) }}
.readiness,.run,.case {{ background:var(--panel); border:1px solid var(--line); border-radius:12px }} .readiness,.run {{ padding:20px; margin:0 0 20px }}
.readiness-grid,.case-grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(320px,1fr)); gap:18px }}
.state {{ padding:10px 12px; border-left:4px solid var(--good); background:#10271f; margin:8px 0 }} .state.missing {{ border-color:var(--warn); background:#2b2114 }}
.case {{ margin:12px 0; overflow:hidden }} summary {{ cursor:pointer; padding:14px; font-weight:700 }} .case>div {{ padding:0 14px 14px }}
.badge {{ font-size:11px; padding:3px 7px; border-radius:999px; margin-right:6px }} .badge.pass {{ color:#052015; background:var(--good) }} .badge.warn {{ color:#291900; background:var(--warn) }}
.metrics {{ border-collapse:collapse; width:100% }} .metrics th,.metrics td {{ text-align:left; padding:5px 8px; border-bottom:1px solid var(--line) }} .metrics th {{ color:var(--muted); font-weight:500 }}
.gallery {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(280px,1fr)); gap:12px; margin-top:14px }} figure {{ margin:0; background:#091018; border-radius:8px; overflow:hidden }} img {{ width:100%; height:240px; object-fit:contain; display:block }} figcaption {{ padding:7px 10px; color:var(--muted) }}
a {{ color:var(--accent) }} .file {{ display:inline-block; padding:3px 7px; margin:2px; border:1px solid var(--line); border-radius:6px; text-decoration:none }} .pass-text {{ color:var(--good) }} .warnings {{ color:var(--warn) }}
</style></head><body>
<header><h1>SHS Physics Results</h1><p class="muted">Generated from {html.escape(str(results_root))}. Rebuild after simulations with <code>python tools/build_results_dashboard.py</code>.</p></header>
<main>
<section class="readiness"><h2>SQUID readiness</h2><div class="readiness-grid">
<div><h3>Validated foundations</h3>
<div class="state">Multiply connected superconducting domains and insulating holes</div>
<div class="state">Gauge-covariant TDGL, fluxoid measurement, current and heat flow around a perforation</div>
<div class="state">Applied field and first self-consistent thin-film magnetic screening implementation</div></div>
<div><h3>Required before a quantitative SQUID claim</h3>
<div class="state missing">Two Josephson weak links or validated interface junction conditions</div>
<div class="state missing">Gauge-invariant phase drop and current–phase relation at each junction</div>
<div class="state missing">Phi0-periodic critical-current modulation and independent reference comparison</div>
<div class="state missing">Circuit bias/readout model, junction capacitance/resistance where required, and noise for switching studies</div></div>
</div></section>
{''.join(run_blocks) if run_blocks else '<section class="run"><p>No run summaries found.</p></section>'}
</main></body></html>"""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(page, encoding="utf-8")
    return len(summaries)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path("benchmark_results"))
    parser.add_argument("--output", type=Path, default=Path("benchmark_results/index.html"))
    args = parser.parse_args()
    count = build_dashboard(args.results, args.output)
    print(f"Wrote {args.output} with {count} result run(s).")


if __name__ == "__main__":
    main()
