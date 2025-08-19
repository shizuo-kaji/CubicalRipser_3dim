"""Python package for CubicalRipser bindings and utilities.

This package exposes the pybind implementation from ``. _cripser`` as a
stable import path ``cripser.computePH`` and provides helpers in
``cripser.utils``.
"""

from .utils import *
from ._cripser import computePH, __version__  # type: ignore
try:
    from tcripser import computePH as computePH_T
except ImportError:
    ValueError(
        "tcripser is not installed. Please install it to use the T-construction."
    )

__all__ = ["computePH", "computePH_T",
    "__version__", "compute_ph",
    "to_gudhi_diagrams",
    "to_gudhi_persistence",
    "group_by_dim"]
