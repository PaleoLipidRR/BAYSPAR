"""Exception types. Messages are MATLAB's, verbatim (PORTING.md GEN-07)."""


class BaysparError(Exception):
    """Base class for every error this package raises."""


class SearchToleranceError(BaysparError, ValueError):
    """No calibration cell falls within the analogue search tolerance."""

    def __init__(self, message: str = "Your search tolerance is too narrow") -> None:
        super().__init__(message)


class EnsembleSizeError(BaysparError, ValueError):
    """More draws were requested than the fitted ensemble holds.

    Args:
        available: Draws in the parameter store.
        requested: Draws the caller asked for.
    """

    def __init__(self, available: int, requested: int) -> None:
        self.available = available
        self.requested = requested
        super().__init__(
            f"requested {requested} draws but the ensemble holds {available}"
        )


class AmbiguousGridCellError(BaysparError, ValueError):
    """A site falls exactly on a calibration cell boundary (PORTING.md BT-08)."""
