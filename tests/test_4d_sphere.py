import math
import numpy as np
import cripser
import signal


def test_4d_sphere_ph():
    # Compute up to H3 on the 4D array via Python module
    arr = np.load('sample/4d_hole.npy').astype(np.float64)
    ph = cripser.compute_ph(arr, maxdim=3, filtration="V")
    check(ph)
    pht = cripser.compute_ph(arr, maxdim=3, filtration="T")
    check(pht)

def check(ph):
    dims = ph[:, 0].astype(int)
    births = ph[:, 1]
    deaths = ph[:, 2]

    # Exactly one H0 infinite bar at birth 0.0
    h0_inf_idx = np.where((dims == 0) & (deaths>1e20))[0]
    assert len(h0_inf_idx) == 1
    assert abs(births[h0_inf_idx[0]] - 0.0) < 1e-9

    # Exactly one H3 bar (0.0, 1.0)
    h3_idx = np.where(dims == 3)[0]
    assert len(h3_idx) == 1
    b = float(births[h3_idx[0]])
    d = float(deaths[h3_idx[0]])
    assert abs(b - 0.0) < 1e-9
    assert abs(d - 1.0) < 1e-9
