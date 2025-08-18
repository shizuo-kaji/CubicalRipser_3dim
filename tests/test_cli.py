import subprocess
import os
import glob
import sys

import pytest


@pytest.fixture(scope="session", autouse=True)
def build_execs():
    pattern = os.path.join('build', 'temp.*', 'cubicalripser')
    if not glob.glob(pattern):
        subprocess.run([sys.executable, 'setup.py', 'build_ext', '--inplace'], check=True)

def get_exec(name):
    pattern = os.path.join('build', 'temp.*', name)
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(f'Executable {name} not found')
    return matches[0]

def test_cubicalripser_cli(tmp_path):
    exe = get_exec('cubicalripser')
    out = tmp_path / 'out.csv'
    subprocess.run([exe, '--output', str(out), 'sample/3dimsample.txt'], check=True)
    with out.open() as f:
        lines = f.readlines()
    assert len(lines) == 18

def test_tcubicalripser_cli(tmp_path):
    exe = get_exec('tcubicalripser')
    out = tmp_path / 'out.csv'
    subprocess.run([exe, '--output', str(out), 'sample/3dimsample.txt'], check=True)
    with out.open() as f:
        lines = f.readlines()
    assert len(lines) == 14
