from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize(["states.pyx", "calculate_error.pyx", "detonations.pyx"])
)
