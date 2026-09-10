import numpy as np
import pytest

from shs.numerics.convergence import (
    ConvergenceController,
    ConvergenceMonitor,
    ConvergenceStatus,
)


def test_monitor_starts_without_history():

    monitor = ConvergenceMonitor()

    assessment = monitor.assess()

    assert assessment.status == ConvergenceStatus.START
    assert assessment.should_continue


def test_monitor_detects_convergence():

    monitor = ConvergenceMonitor(
        tolerance=1e-6,
    )

    assessment = monitor.record(
        iteration=1,
        residual=1e-5,
    )

    assert assessment.status == ConvergenceStatus.START

    assessment = monitor.record(
        iteration=2,
        residual=1e-7,
    )

    assert assessment.status == ConvergenceStatus.CONVERGED
    assert assessment.should_continue is False
    assert assessment.estimated_iterations_remaining == 0.0


def test_monitor_detects_fast_convergence():

    monitor = ConvergenceMonitor(
        tolerance=1e-12,
    )

    monitor.record(1, 1e-2)
    assessment = monitor.record(2, 1e-3)

    assert assessment.status == ConvergenceStatus.FAST


def test_monitor_detects_healthy_convergence():

    monitor = ConvergenceMonitor(
        tolerance=1e-12,
    )

    monitor.record(1, 1e-2)
    assessment = monitor.record(2, 7e-3)

    assert assessment.status == ConvergenceStatus.HEALTHY


def test_monitor_detects_slow_convergence():

    monitor = ConvergenceMonitor(
        tolerance=1e-12,
    )

    monitor.record(1, 1e-2)
    assessment = monitor.record(2, 9.8e-3)

    assert assessment.status == ConvergenceStatus.SLOW


def test_monitor_detects_stalling():

    monitor = ConvergenceMonitor(
        tolerance=1e-12,
    )

    monitor.record(1, 1e-3)
    assessment = monitor.record(2, 1e-3)

    assert assessment.status == ConvergenceStatus.STALLING
    assert assessment.adaptation_recommended


def test_monitor_detects_divergence():

    monitor = ConvergenceMonitor(
        tolerance=1e-12,
    )

    monitor.record(1, 1e-3)
    assessment = monitor.record(2, 2e-3)

    assert assessment.status == ConvergenceStatus.DIVERGING
    assert assessment.adaptation_recommended


def test_monitor_detects_oscillation():

    monitor = ConvergenceMonitor(
        tolerance=1e-12,
    )

    monitor.record(1, 1.0)
    monitor.record(2, 0.5)
    monitor.record(3, 0.6)

    assessment = monitor.record(
        4,
        0.3,
    )

    assert assessment.status == (
        ConvergenceStatus.OSCILLATING
    )
    assert assessment.adaptation_recommended


def test_iteration_prediction():

    monitor = ConvergenceMonitor(
        tolerance=1e-8,
    )

    for iteration in range(1, 11):
        residual = 10.0 ** (-iteration / 2.0)
        monitor.record(
            iteration,
            residual,
        )

    estimate = monitor.estimate_iterations_remaining()

    assert estimate is not None
    assert estimate > 0.0


def test_prediction_confidence_is_bounded():

    monitor = ConvergenceMonitor()

    for iteration in range(1, 8):
        monitor.record(
            iteration,
            10.0 ** (-iteration / 2.0),
        )

    confidence = monitor.prediction_confidence()

    assert 0.0 <= confidence <= 1.0


def test_checkpoint_interval():

    monitor = ConvergenceMonitor()

    controller = ConvergenceController(
        monitor=monitor,
        check_interval=10,
    )

    assert controller.should_check(1)
    assert not controller.should_check(5)
    assert controller.should_check(10)
    assert controller.should_check(20)


def test_controller_enforces_max_iterations():

    monitor = ConvergenceMonitor(
        tolerance=1e-20,
    )

    controller = ConvergenceController(
        monitor=monitor,
        max_iterations=5,
    )

    controller.record(1, 1e-2)
    controller.record(2, 1e-2)
    controller.record(3, 1e-2)
    controller.record(4, 1e-2)
    controller.record(5, 1e-2)

    assessment = controller.decision(5)

    assert assessment.status == (
        ConvergenceStatus.MAX_ITERATIONS
    )
    assert assessment.should_continue is False
    assert assessment.adaptation_recommended


def test_monitor_rejects_invalid_residual():

    monitor = ConvergenceMonitor()

    with pytest.raises(ValueError):
        monitor.record(
            iteration=1,
            residual=-1.0,
        )


def test_monitor_rejects_nonfinite_residual():

    monitor = ConvergenceMonitor()

    with pytest.raises(ValueError):
        monitor.record(
            iteration=1,
            residual=np.inf,
        )


def test_monitor_requires_increasing_iterations():

    monitor = ConvergenceMonitor()

    monitor.record(1, 1e-3)

    with pytest.raises(ValueError):
        monitor.record(1, 1e-4)