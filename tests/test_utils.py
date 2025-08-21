import numpy as np

from cripser import (
    compute_ph,
    to_gudhi_diagrams,
    to_gudhi_persistence,
    group_by_dim,
)


def test_utils_converters_basic():
    # Simple 2D constant image ensures an infinite H0 bar
    arr = np.zeros((3, 3), dtype=np.float64)
    ph = compute_ph(arr, maxdim=1)

    # Base shape checks
    assert ph.ndim == 2 and ph.shape[1] == 9

    # Diagrams
    dgms = to_gudhi_diagrams(ph, maxdim=1)
    assert isinstance(dgms, list) and len(dgms) == 2
    assert dgms[0].ndim == 2 and dgms[0].shape[1] == 2
    # Expect at least one infinite death in H0
    assert np.isinf(dgms[0][:, 1]).any()

    # Persistence list
    pers = to_gudhi_persistence(ph)
    assert any(d == 0 and np.isinf(bd[1]) for d, bd in pers)

    # Group by dimension
    groups = group_by_dim(ph)
    assert len(groups) >= 1
    if len(groups[0]) > 0:
        assert np.all(groups[0][:, 0] == 0)
