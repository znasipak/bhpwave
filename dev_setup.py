import os
from Cython.Build import cythonize
from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np
import sys


"""
For some reason, the libraries extension does not support extra compiler arguments. So
instead, we pass them directly by changing the evironment variable CFLAGS
"""
compiler_flags = ["-std=c++11", "-fopenmp", "-march=native"]
libraries = ["gsl", "gslcblas"]

# set some flags and libraries depending on the system platform
if sys.platform.startswith('win32'):
    compiler_flags.append('/Od')
    libraries.append('gomp')
elif sys.platform.startswith('darwin'):
    compiler_flags.append('-O2')
    libraries.append('omp')
elif sys.platform.startswith('linux'):
    compiler_flags.append('-O2')
    libraries.append('gomp')
    

CFLAGS = os.getenv("CFLAGS")
if CFLAGS is None:
    CFLAGS = ""
else:
    CFLAGS = str(CFLAGS)
os.environ["CFLAGS"] = CFLAGS
os.environ["CFLAGS"] += " "
os.environ["CFLAGS"] += ' '.join(compiler_flags)
base_path = sys.prefix

spline_dependence = ["cpp/src/spline.cpp"]
swsh_dependence = ["cpp/src/swsh.cpp"]
trajectory_dependence = ["cpp/src/trajectory.cpp", *spline_dependence]
harmonics_dependence = ["cpp/src/harmonics.cpp", *spline_dependence]
waveform_dependence = ["cpp/src/waveform.cpp", *harmonics_dependence, *trajectory_dependence, *swsh_dependence]
fourier_dependence = ["cpp/src/fourier.cpp", *waveform_dependence]

full_dependence = [*fourier_dependence]

lib_extension = dict(
    sources = [*set(full_dependence)],
    # libraries=["gsl", "gslcblas", "lapack", "lapacke", "omp"],
    libraries=libraries,
    language='c++',
    include_dirs = ["cpp/include", base_path + "/include"],
)
bhpwavecpp = ['bhpwavecpp', lib_extension]

cpu_extension = dict(
    libraries=libraries,
    language='c++',
    include_dirs=["cpp/include", np.get_include(), base_path + "/include"],
)

swsh_ext = Extension(
    "bhpswshcy", 
    sources=["cython/swsh_wrap.pyx"], 
    **cpu_extension,
)

wave_ext = Extension(
    "bhpwaveformcy", 
    sources=["cython/waveform_wrap.pyx"], 
    **cpu_extension,
)

ext_modules = [swsh_ext, wave_ext]

setup(
    name = "bhpwave",
    author = "Zach Nasipak",
    version = "0.1.0",
    description = "Adiabatic EMRI waveform generator",
    ext_modules = cythonize(ext_modules, language_level = "3"),
    packages = ["bhpwave", "bhpwave.swsh", "bhpwave.trajectory", "bhpwave.harmonics", "bhpwave.data"],
    py_modules = ["bhpwave.waveform", "bhpwave.constants", "bhpwave.spline"],
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Natural Language :: English",
        "Programming Language :: C++",
        "Programming Language :: Cython",
    ],
    libraries = [bhpwavecpp],
    cmdclass = {'build_ext': build_ext},
    package_data = {"bhpwave.data": ["*.txt"]},
    include_package_data = True,
    zip_safe = False
)