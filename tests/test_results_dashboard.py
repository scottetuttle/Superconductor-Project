from pathlib import Path

from tools.build_results_dashboard import build_dashboard


def test_dashboard_embeds_metrics_media_and_sqid_readiness():
    results = Path("benchmark_results")
    output = results / "index.html"
    assert build_dashboard(results, output) >= 1
    page = output.read_text(encoding="utf-8")
    assert "SQUID readiness" in page
    assert "transport current A" in page
    assert "diagnostics.csv" in page
