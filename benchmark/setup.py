# setup.py — 编译Cython代码
from setuptools import setup
from Cython.Build import cythonize
import numpy as np

setup(
    ext_modules=cythonize(
        "normalize_cython.pyx",
        compiler_directives={
            'language_level': '3',
            'boundscheck': False,
            'wraparound': False,
        }
    ),
    include_dirs=[np.get_include()]
)
