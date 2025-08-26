import numpy as np

import cripser
import tcripser
import gudhi as gd

def _against_gudhi(dim=4,filtration="V"):
    arr = np.random.random([10]*dim)
    ph = cripser.compute_ph(arr, maxdim=dim-1, filtration=filtration)
    ph = cripser.to_gudhi_diagrams(ph)
    if filtration=="V":
        cc = gd.CubicalComplex(vertices=arr)
    else:
        cc = gd.CubicalComplex(top_dimensional_cells=arr)
    cc.persistence(homology_coeff_field=2, min_persistence=0)
    for d in range(dim-1):
        cripser_ints =  np.asarray(ph[d])
        gudhi_ints = np.asarray(cc.persistence_intervals_in_dimension(d))
        # Sort intervals lexicographically by (birth, death) for stable comparison
        def sort_ints(a):
            if a.size == 0:
                return a.reshape(0, 2)
            return a[np.lexsort((a[:, 1], a[:, 0]))]
        ci = sort_ints(cripser_ints)
        gi = sort_ints(gudhi_ints)
        print("cripser", d,ci.shape)
        print(gi.shape)
        assert ci.shape == gi.shape, f"Mismatch count in dim {d}: {ci.shape} vs {gi.shape}"
        assert np.allclose(ci, gi), f"Intervals differ in dim {d}"

def test_cripser_vs_gudhi():
    for d in range(1, 4):
        for filtration in ["V", "T"]:
           _against_gudhi(dim=d, filtration=filtration)
