"""
Convergence monitoring and control infrastructure for SHS.

This module provides numerical tools for monitoring iterative
convergence independently of any particular physics model.

The monitor observes residual histories and classifies their
behavior. The controller interprets that information and decides
whether an iterative process should continue, accept convergence,
or request numerical adaptation.

The controller is intentionally observational at this stage.
It does not modify solver parameters automatically.
"""

from dataclasses import dataclass
from enum import Enum

import numpy as np


class ConvergenceStatus(Enum):
    """
    Classification of the current convergence behavior.
    """

    START = "start"
    FAST = "fast"
    HEALTHY = "healthy"
    SLOW = "slow"
    STALLING = "stalling"
    OSCILLATING = "oscillating"
    DIVERGING = "diverging"
    CONVERGED = "converged"
    MAX_ITERATIONS = "max_iterations"


@dataclass
class ConvergenceAssessment:
    """
    Snapshot of the current convergence behavior.
    """

    iteration: int
    residual: float
    tolerance: float
    status: ConvergenceStatus

    residual_ratio: float | None = None
    log_slope: float | None = None

    estimated_iterations_remaining: float | None = None
    prediction_confidence: float = 0.0

    should_continue: bool = True
    adaptation_recommended: bool = False


@dataclass
class ConvergenceMonitor:
    """
    Monitors residual history for an iterative numerical process.

    The monitor does not modify the process being monitored.
    """

    tolerance: float = 1e-8
    prediction_window: int = 10

    fast_ratio: float = 0.5
    healthy_ratio: float = 0.95
    stall_ratio: float = 0.999

    reasonable_iterations: int = 50

    def __post_init__(self):
        self.validate()

        self.iterations = []
        self.residuals = []

    def validate(self):
        """
        Validate convergence-monitor configuration.
        """

        if not np.isfinite(self.tolerance):
            raise ValueError(
                "Convergence tolerance must be finite."
            )

        if self.tolerance <= 0.0:
            raise ValueError(
                "Convergence tolerance must be positive."
            )

        if self.prediction_window < 2:
            raise ValueError(
                "Prediction window must be at least 2."
            )

        if not (
            0.0 < self.fast_ratio < 1.0
        ):
            raise ValueError(
                "Fast ratio must lie between 0 and 1."
            )

        if not (
            0.0 < self.healthy_ratio < 1.0
        ):
            raise ValueError(
                "Healthy ratio must lie between 0 and 1."
            )

        if not (
            0.0 < self.stall_ratio <= 1.0
        ):
            raise ValueError(
                "Stall ratio must lie between 0 and 1."
            )

        if self.reasonable_iterations < 1:
            raise ValueError(
                "Reasonable iterations must be positive."
            )

        return True

    def reset(self):
        """
        Clear all stored convergence history.
        """

        self.iterations.clear()
        self.residuals.clear()

    @property
    def iteration(self):
        """
        Return the most recent iteration number.
        """

        if not self.iterations:
            return 0

        return self.iterations[-1]

    @property
    def residual(self):
        """
        Return the most recent residual.
        """

        if not self.residuals:
            return np.inf

        return self.residuals[-1]

    def record(
        self,
        iteration: int,
        residual: float,
    ):
        """
        Record one residual observation.

        Parameters
        ----------
        iteration:
            Iteration number.

        residual:
            Non-negative convergence residual.
        """

        iteration = int(iteration)
        residual = float(residual)

        if iteration < 0:
            raise ValueError(
                "Iteration number cannot be negative."
            )

        if not np.isfinite(residual):
            raise ValueError(
                "Residual must be finite."
            )

        if residual < 0.0:
            raise ValueError(
                "Residual cannot be negative."
            )

        if self.iterations:
            if iteration <= self.iterations[-1]:
                raise ValueError(
                    "Iteration numbers must increase."
                )

        self.iterations.append(iteration)
        self.residuals.append(residual)

        return self.assess()

    def residual_ratio(self):
        """
        Return the ratio of the latest residual to the previous
        residual.

        A value below one indicates decreasing residual.
        """

        if len(self.residuals) < 2:
            return None

        previous = self.residuals[-2]

        if previous == 0.0:
            return 0.0

        return self.residuals[-1] / previous

    def log_slope(self):
        """
        Estimate the slope of log10(residual) over the most recent
        prediction window.

        Negative values indicate convergence.
        """

        if len(self.residuals) < 2:
            return None

        count = min(
            len(self.residuals),
            self.prediction_window,
        )

        residuals = np.asarray(
            self.residuals[-count:],
            dtype=float,
        )

        iterations = np.asarray(
            self.iterations[-count:],
            dtype=float,
        )

        if np.any(residuals <= 0.0):
            return None

        log_residuals = np.log10(residuals)

        slope, _ = np.polyfit(
            iterations,
            log_residuals,
            1,
        )

        return float(slope)

    def estimate_iterations_remaining(self):
        """
        Estimate how many additional iterations are required
        to reach the convergence tolerance.

        The estimate is based on the recent linear trend in
        log10(residual).
        """

        residual = self.residual

        if residual <= self.tolerance:
            return 0.0

        slope = self.log_slope()

        if slope is None:
            return None

        if slope >= 0.0:
            return None

        log_distance = (
            np.log10(self.tolerance)
            - np.log10(residual)
        )

        estimate = log_distance / slope

        if not np.isfinite(estimate):
            return None

        estimate = max(float(estimate), 0.0)

        return estimate

    def prediction_confidence(self):
        """
        Estimate confidence in the iteration prediction.

        Confidence is based on the amount of history available
        and the consistency of the recent convergence trend.
        """

        if len(self.residuals) < 2:
            return 0.0

        count = min(
            len(self.residuals),
            self.prediction_window,
        )

        if count < 3:
            return 0.25

        residuals = np.asarray(
            self.residuals[-count:],
            dtype=float,
        )

        if np.any(residuals <= 0.0):
            return 0.0

        iterations = np.asarray(
            self.iterations[-count:],
            dtype=float,
        )

        log_residuals = np.log10(residuals)

        coefficients = np.polyfit(
            iterations,
            log_residuals,
            1,
        )

        predicted = np.polyval(
            coefficients,
            iterations,
        )

        actual = log_residuals

        variance = np.sum(
            (actual - np.mean(actual)) ** 2
        )

        if variance == 0.0:
            return 1.0

        error = np.sum(
            (actual - predicted) ** 2
        )

        r_squared = 1.0 - error / variance

        return float(
            np.clip(r_squared, 0.0, 1.0)
        )

    def _is_oscillating(self):
        """
        Detect alternating increases and decreases in residual.
        """

        if len(self.residuals) < 4:
            return False

        recent = np.asarray(
            self.residuals[-4:],
            dtype=float,
        )

        differences = np.diff(recent)

        signs = np.sign(differences)

        signs = signs[signs != 0.0]

        if len(signs) < 3:
            return False

        alternating = np.all(
            signs[1:] != signs[:-1]
        )

        return bool(alternating)

    def assess(self):
        """
        Assess the current convergence behavior.
        """

        iteration = self.iteration
        residual = self.residual

        if (
            residual <= self.tolerance
            and iteration >= 5
        ):

            return ConvergenceAssessment(
                iteration=iteration,
                residual=residual,
                tolerance=self.tolerance,
                status=ConvergenceStatus.CONVERGED,
                residual_ratio=self.residual_ratio(),
                log_slope=self.log_slope(),
                estimated_iterations_remaining=0,
                prediction_confidence=1.0,
                should_continue=False,
                adaptation_recommended=False,
            )

        if len(self.residuals) < 2:
            return ConvergenceAssessment(
                iteration=iteration,
                residual=residual,
                tolerance=self.tolerance,
                status=ConvergenceStatus.START,
                should_continue=True,
                adaptation_recommended=False,
            )

        ratio = self.residual_ratio()
        slope = self.log_slope()

        if self._is_oscillating():
            status = ConvergenceStatus.OSCILLATING

        elif ratio is not None and ratio > 1.0:
            status = ConvergenceStatus.DIVERGING

        elif (
            ratio is not None
            and ratio >= self.stall_ratio
        ):
            status = ConvergenceStatus.STALLING

        elif (
            ratio is not None
            and ratio <= self.fast_ratio
        ):
            status = ConvergenceStatus.FAST

        elif (
            ratio is not None
            and ratio <= self.healthy_ratio
        ):
            status = ConvergenceStatus.HEALTHY

        else:
            status = ConvergenceStatus.SLOW

        estimate = self.estimate_iterations_remaining()

        if estimate is not None:
            estimate = max(5.0, estimate)

        confidence = self.prediction_confidence()

        adaptation_recommended = status in {
            ConvergenceStatus.STALLING,
            ConvergenceStatus.OSCILLATING,
            ConvergenceStatus.DIVERGING,
        }

        if (
            estimate is not None
            and estimate > self.reasonable_iterations
        ):
            adaptation_recommended = True

        return ConvergenceAssessment(
            iteration=iteration,
            residual=residual,
            tolerance=self.tolerance,
            status=status,
            residual_ratio=ratio,
            log_slope=slope,
            estimated_iterations_remaining=estimate,
            prediction_confidence=confidence,
            should_continue=True,
            adaptation_recommended=adaptation_recommended,
        )


@dataclass
class ConvergenceController:
    """
    Controls the decision-making layer around a convergence monitor.

    At this stage the controller is intentionally observational.
    It does not modify numerical parameters.
    """

    monitor: ConvergenceMonitor
    max_iterations: int = 1000
    check_interval: int = 10

    def __post_init__(self):
        self.validate()

    def validate(self):
        """
        Validate controller configuration.
        """

        if self.max_iterations < 1:
            raise ValueError(
                "Maximum iterations must be positive."
            )

        if self.check_interval < 1:
            raise ValueError(
                "Check interval must be positive."
            )

        return True

    def reset(self):
        """
        Reset the underlying convergence monitor.
        """

        self.monitor.reset()

    def record(
        self,
        iteration: int,
        residual: float,
    ):
        """
        Record a residual and return the current assessment.
        """

        return self.monitor.record(
            iteration=iteration,
            residual=residual,
        )

    def should_check(self, iteration: int):
        """
        Return whether this iteration is a strategic checkpoint.

        Convergence itself is still assessed every iteration.
        """

        iteration = int(iteration)

        return (
            iteration == 1
            or iteration % self.check_interval == 0
        )

    def decision(self, iteration: int):
        """
        Return the current convergence assessment.

        This method does not alter numerical parameters.
        """

        assessment = self.monitor.assess()

        if assessment.status == (
            ConvergenceStatus.CONVERGED
        ):
            return assessment

        if iteration >= self.max_iterations:
            assessment.status = (
                ConvergenceStatus.MAX_ITERATIONS
            )
            assessment.should_continue = False
            assessment.adaptation_recommended = True

        return assessment