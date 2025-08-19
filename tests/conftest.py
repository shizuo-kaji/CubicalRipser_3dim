import glob
import subprocess
import sys

import pytest


@pytest.fixture(scope="session", autouse=True)
def build_python_modules():
    built = glob.glob('cripser/_cripser*.so') and glob.glob('tcripser*.so')
    if not built:
        subprocess.run([sys.executable, 'setup.py', 'build_ext', '--inplace'], check=True)
