import numpy as np

import cripser
import tcripser


def test_cripser_module_on_3d_hole():
    arr = np.load('sample/3d_hole.npy').astype(np.float64)
    ph = cripser.computePH(arr, maxdim=2)
    assert ph.ndim == 2 and ph.shape[1] == 9
    dims = set(ph[:, 0].astype(int))
    # Expect H0 and H2 features on a 3D hole dataset
    assert 0 in dims
    assert 2 in dims


def test_tcripser_module_on_3d_hole():
    arr = np.load('sample/3d_hole.npy').astype(np.float64)
    ph = tcripser.computePH(arr, maxdim=2)
    assert ph.ndim == 2 and ph.shape[1] == 9
    dims = set(ph[:, 0].astype(int))
    # Expect non-empty result and some higher-dimensional features
    assert len(dims) >= 1
    assert any(d in dims for d in (1, 2))
