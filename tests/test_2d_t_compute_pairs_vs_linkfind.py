import os
import subprocess
import sys
from pathlib import Path

import numpy as np


def run_tcubicalripser(args, cwd=None):
    exe = Path('build') / 'tcubicalripser'
    if not exe.exists():
        # Try to build the binary if not present
        subprocess.run(['make', '-C', 'build', 'tcubicalripser'], check=True)
    cmd = [str(exe)] + args
    proc = subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)
    return proc.stdout


def load_csv(path):
    # CSV columns: dim,birth,death,(coords...)
    data = np.loadtxt(path, delimiter=',')
    if data.ndim == 1:
        data = data[None, :]
    return data


def test_2d_t_compute_pairs_matches_linkfind(tmp_path):
    # Deterministic 2D test image (16x16): smooth gradient with a bump
    x = np.linspace(0, 1, 16)
    y = np.linspace(0, 1, 16)
    X, Y = np.meshgrid(x, y, indexing='xy')
    img = (X + Y) / 2.0
    # Add a small localized feature to create nontrivial PH1
    cx, cy = 8, 8
    rr = (np.arange(16) - cx)[None, :] ** 2 + (np.arange(16)[:, None] - cy) ** 2
    img += 0.05 * np.exp(-rr / 8.0)
    img = img.astype(np.float64)

    # random
    img = np.random.RandomState(0).rand(16, 16)
    # Save as .npy (tcubicalripser supports .npy input)
    npy_path = tmp_path / 'img.npy'
    np.save(npy_path, img)

    out_link = tmp_path / 'link.csv'
    out_comp = tmp_path / 'comp.csv'

    # Limit to maxdim=1 for this test
    run_tcubicalripser(['-m', '1', '-o', str(out_link), str(npy_path)])
    run_tcubicalripser(['-m', '1', '-a', 'compute_pairs', '-o', str(out_comp), str(npy_path)])

    A = load_csv(out_link)
    B = load_csv(out_comp)

    # Assert same number of pairs in H0 and H1
    for d in (0, 1):
        a = A[A[:, 0] == d]
        b = B[B[:, 0] == d]
        assert a.shape[0] == b.shape[0], f"dim {d} count mismatch: {a.shape[0]} vs {b.shape[0]}"

    # Compare full rows up to columns present (allow exact float equality for deterministic data)
    # Sort by (dim, birth, death, x1,y1,x2,y2)
    def sort_key(rows):
        return rows[np.lexsort(np.flipud(rows[:, :3].T))]

    A_sorted = A[np.lexsort((A[:, 2], A[:, 1], A[:, 0]))]
    B_sorted = B[np.lexsort((B[:, 2], B[:, 1], B[:, 0]))]

    assert A_sorted.shape == B_sorted.shape

    # Direct comparison; if this is too strict on some platforms, consider np.allclose with tight tol
    if not np.array_equal(A_sorted, B_sorted):
        # Provide a brief diff summary for debugging
        diff_count = np.sum(~np.isclose(A_sorted, B_sorted))
        raise AssertionError(f"CSV mismatch: {diff_count} differing entries\nA[:5]={A_sorted[:5]}\nB[:5]={B_sorted[:5]}")
