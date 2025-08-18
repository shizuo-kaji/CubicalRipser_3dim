import csv
import glob
import math
import os
import subprocess


def get_exec(name: str) -> str:
    pattern = os.path.join('build', 'temp.*', name)
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(f'Executable {name} not found')
    return matches[0]


def read_csv_rows(path: str):
    rows = []
    with open(path, newline='') as f:
        r = csv.reader(f)
        for row in r:
            if not row:
                continue
            # convert first three columns: dim, birth, death
            dim = int(float(row[0]))
            birth = float(row[1])
            death = float(row[2])
            rows.append((dim, birth, death))
    return rows


def test_4d_sphere_ph(tmp_path):
    exe = get_exec('cubicalripser')
    out = tmp_path / 'out.csv'

    # Compute up to H3 on the 4D array
    subprocess.run([exe, '--maxdim', '3', '--output', str(out), 'sample/4d_hole.npy'], check=True)

    rows = read_csv_rows(str(out))

    # Expect exactly one H0 infinite bar at birth 0.0
    h0_inf = [(b, d) for (dim, b, d) in rows if dim == 0 and math.isinf(d)]
    assert len(h0_inf) == 1
    assert abs(h0_inf[0][0] - 0.0) < 1e-9

    # Expect exactly one H3 bar (0.0, 1.0)
    h3 = [(b, d) for (dim, b, d) in rows if dim == 3]
    assert len(h3) == 1
    b, d = h3[0]
    assert abs(b - 0.0) < 1e-9
    assert abs(d - 1.0) < 1e-9
