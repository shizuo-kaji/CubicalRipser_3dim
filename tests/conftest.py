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

