from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EDITOR = ROOT / "tools" / "config_editor"


def test_config_editor_assets_and_export_contract_exist():
    html = (EDITOR / "index.html").read_text(encoding="utf-8")
    javascript = (EDITOR / "app.js").read_text(encoding="utf-8")
    assert "deviceCanvas" in html
    assert "Circle hole" in html
    assert "+ Vortex" in html
    assert "Laser path" in html
    assert "Optical tweezer" in html
    assert "SNS / weak-link junction" in html
    assert "junctionControls" in html
    assert "feedbackControls" in html
    for function in ("geometry()", "simulation()", "diagnostics()", "project()"):
        assert function in javascript
    assert "const PRESETS" in javascript
    assert "optical_junction" in javascript
    assert "applyPreset()" in javascript
    assert 'runner:"feedback"' in javascript
    assert "run_command:" in javascript
    assert "_experiment.json" in javascript
    assert 'power_fraction:p.power' in javascript
    assert 'x_fraction:v.x/n("widthNm")' in javascript
    assert 'physics_backend' not in javascript  # Exports the full TDGL path.


def test_editor_server_uses_loopback_by_default():
    server = (ROOT / "tools" / "serve_config_editor.py").read_text(encoding="utf-8")
    assert 'default="127.0.0.1"' in server
    assert 'default=8765' in server
