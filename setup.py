import setuptools  # noqa: F401
from numpy.distutils.core import setup, Extension

EXT1 = Extension(name='propagate', sources=['propagate.f90'],
                 extra_f90_compile_args=['-fopenmp'],
                 extra_compile_args=['-lgomp'])

setup(ext_modules=[EXT1])
