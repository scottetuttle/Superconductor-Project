"""Low-overhead terminal progress reporting for step-based simulations."""

from __future__ import annotations

import math
import time


class PercentProgress:
    """Print each newly reached integer percentage and an estimated finish time."""

    def __init__(self, total_steps, *, label="run", completed_steps=0, stream=None,
                 clock=time.monotonic):
        if not isinstance(total_steps, int) or total_steps < 1:
            raise ValueError("total_steps must be a positive integer.")
        if completed_steps < 0 or completed_steps > total_steps:
            raise ValueError("completed_steps must be between zero and total_steps.")
        self.total_steps = total_steps
        self.label = str(label)
        self.stream = stream
        self.clock = clock
        self.started = clock()
        self.started_step = int(completed_steps)
        self.last_percent = math.floor(100 * completed_steps / total_steps)

    def update(self, completed_steps):
        """Report all integer percentages crossed by ``completed_steps``."""
        completed = min(max(int(completed_steps), 0), self.total_steps)
        reached = math.floor(100 * completed / self.total_steps)
        if reached <= self.last_percent:
            return
        elapsed = max(self.clock() - self.started, 0.0)
        work = completed - self.started_step
        rate = elapsed / work if work > 0 else 0.0
        eta = rate * (self.total_steps - completed)
        for percent in range(self.last_percent + 1, reached + 1):
            threshold = math.ceil(percent * self.total_steps / 100)
            shown_step = min(max(threshold, self.started_step), completed)
            message = (f"[{self.label}] {percent:3d}% "
                       f"({shown_step}/{self.total_steps} steps) "
                       f"elapsed={elapsed:.1f}s eta={eta:.1f}s")
            print(message, file=self.stream, flush=True)
        self.last_percent = reached
