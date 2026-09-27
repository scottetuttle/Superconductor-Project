"""Safe output-directory allocation for repeatable simulation runners."""
from pathlib import Path


def reserve_output_directory(base):
    """Create base or the first available numbered sibling without overwriting."""
    base = Path(base)
    for index in range(10000):
        candidate = base if index == 0 else base.parent / f"{base.name}_{index}"
        try:
            candidate.mkdir(parents=True, exist_ok=False)
            return candidate
        except FileExistsError:
            continue
    raise RuntimeError(f"Could not reserve a unique output directory for {base}.")
