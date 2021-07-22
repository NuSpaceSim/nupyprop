import setuptools  # noqa: F401
from numpy.distutils.core import setup, Extension

EXT1 = Extension(
    name="propagate",
    sources=["src/nupyprop/propagate.f90"],
    extra_f90_compile_args=["-fopenmp"],
    extra_link_args=["-fopenmp"],
)

setup(ext_modules=[EXT1])
