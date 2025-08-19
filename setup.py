# setup.py
import os, sys, re, platform, subprocess
from distutils.version import LooseVersion
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    pass

if __name__ == "__main__":
    setup(
        ext_modules=[CMakeExtension("cripser"), CMakeExtension("tcripser")],
        cmdclass={"build_ext": CMakeBuild},
    )
