import io

from shs.utils.progress import PercentProgress


def test_reports_each_crossed_integer_percentage_once():
    stream = io.StringIO()
    reporter = PercentProgress(
        250, label="case", stream=stream, clock=lambda: 2.0)

    reporter.update(5)
    reporter.update(5)
    reporter.update(250)

    lines = stream.getvalue().splitlines()
    assert len(lines) == 100
    assert lines[0].startswith("[case]   1%")
    assert lines[1].startswith("[case]   2%")
    assert lines[-1].startswith("[case] 100%")


def test_resume_starts_at_next_unreported_percentage():
    stream = io.StringIO()
    reporter = PercentProgress(
        1000, completed_steps=405, stream=stream, clock=lambda: 1.0)

    reporter.update(410)

    assert stream.getvalue().splitlines()[0].startswith("[run]  41%")
