import os
from Cython.Build import cythonize
from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np
from subprocess import run

compiler_flags = ["-std=c++11", "-fopenmp", "-O2", "-march=native"]

try:
    CC_string = os.getenv("CC")
except:
    os.environ["CC"] = ""

os.environ["CC"] += " "
os.environ["CC"] += ' '.join(compiler_flags)

print(os.environ["CC"])

spline_dependence = ["cpp/src/spline.cpp"]
swsh_dependence = ["cpp/src/swsh.cpp"]
trajectory_dependence = ["cpp/src/trajectory.cpp", *spline_dependence]
harmonics_dependence = ["cpp/src/harmonics.cpp", *spline_dependence]
waveform_dependence = ["cpp/src/waveform.cpp", *harmonics_dependence, *trajectory_dependence, *swsh_dependence]

full_dependence = [*waveform_dependence]

lib_extension = dict(
    sources = [*set(full_dependence)],
    libraries=["gsl", "gslcblas", "lapack", "lapacke", "omp"],
    language='c++',
    include_dirs = ["cpp/include"],
)
bhpwavecpp = ['bhpwavecpp', lib_extension]

cpu_extension = dict(
    libraries=["gsl", "gslcblas", "lapack", "lapacke", "omp"],
    language='c++',
    include_dirs=["cpp/include", np.get_include()],
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
    py_modules = ["bhpwave.waveform"],
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