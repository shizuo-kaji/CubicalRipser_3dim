import os
import re
import sys
import platform
import subprocess
from pathlib import Path

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


class CMakeExtension(Extension):
    """Extension for building C++ code with CMake."""
    def __init__(self, name, sourcedir=''):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    """Custom build class for handling CMake-based builds."""
    
    def run(self):
        # Check if CMake is installed
        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions)
            )

        # Check for minimum CMake version on Windows
        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r'version\s*([\d.]+)', 
                subprocess.check_output(['cmake', '--version']).decode()
            ).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        # Build each extension
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        # Directory for output libraries
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # Ensure the output directory ends with a separator
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        # CMake configuration arguments
        cmake_args = [
            f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}',
            f'-DPYTHON_EXECUTABLE={sys.executable}'
        ]

        # Build configuration
        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        # Platform-specific adjustments
        if platform.system() == "Windows":
            cmake_args += [f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}']
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += [f'-DCMAKE_BUILD_TYPE={cfg}']
            build_args += ['--', '-j2']

        # Prepare environment variables
        env = os.environ.copy()
        env['CXXFLAGS'] = f'{env.get("CXXFLAGS", "")} -DVERSION_INFO=\\"{self.distribution.get_version()}\\"'

        # Ensure build directory exists
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        # Run CMake to configure and build the extension
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)


# Read the README file for long description
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Setup configuration
setup(
    name='cripser',
    version='0.0.14',
    author='Shizuo KAJI',
    author_email='shizuo.kaji@gmail.com',
    description='Cubical Ripser Python binding',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='MIT',
    url='https://github.com/shizuo-kaji/CubicalRipser_3dim',
    keywords='persistent homology TDA topological image volume',
    ext_modules=[
        CMakeExtension('cripser'),
        CMakeExtension('tcripser')
    ],
    cmdclass={'build_ext': CMakeBuild},
    zip_safe=False,
)