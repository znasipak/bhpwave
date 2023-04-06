# bhpwave
An adiabatic waveform generator using black hole perturbation theory

BHPWave generates gravitational waveforms of extreme-mass-ratio inspirals (EMRIs) using the adiabatic approximation of black hole perturbation theory. The model is restricted to binaries in which the small body is undergoing a quasi-circular, equatorial inspiral into a rotating massive black hole.

# Installation

BHPWave relies on a few dependencies to install and run, namely
a C/C++ compiler (e.g., `g++`), `gsl`, `Cython`, 
`numpy`, and `python >= 3.7`.
To reduce package conflicts and ensure that the proper dependencies are installed,
we recommend using Anaconda and its virtual environments.

Create a conda environment `bhpwave-env` (or whatever name you would like)
with the necessary dependencies to install `bhpwave`:
```
$ conda create -n bhpwave-env -c conda-forge gsl Cython numpy python=3.7
$ conda activate bhpwave-env
```
To include the necessary compiler on MACOSX run:
```
$ conda install clang_osx-64 clangxx_osx-64
```
To include the necessary compiler on linux run:
```
$ conda install gcc_linux-64 gxx_linux-64
```
Next clone the :code:`bhpwave` repository from GitHub:
```
$ git clone https://github.com/znasipak/bhpwave.git
$ cd bhpwave
```
Finally, we recommend installing the package via `pip`:
```
$ pip install .
```
