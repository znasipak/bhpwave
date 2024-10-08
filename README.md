# bhpwave
[![Documentation Status](https://readthedocs.org/projects/bhpwave/badge/?version=latest)](https://bhpwave.readthedocs.io/en/latest/?badge=latest)
[![Read the article](https://img.shields.io/badge/PDF-latest-blue.svg?style=flat)](https://github.com/znasipak/bhpwave-article/raw/gh-action-result/pdflatex/ms.pdf)
[![DOI](https://img.shields.io/badge/arXiv-2310.19706-B31B1B)](https://doi.org/10.48550/arXiv.2310.19706)
[![DOI](https://img.shields.io/badge/PhysRevD-109.044020-purple)](https://doi.org/10.1103/PhysRevD.109.044020)
[![DOI](https://img.shields.io/badge/znasipak-bhpwave--article-blue?logo=Github)](https://github.com/znasipak/bhpwave-article)
[![DOI](https://zenodo.org/badge/621884293.svg)](https://zenodo.org/doi/10.5281/zenodo.13289544)



An adiabatic waveform generator using black hole perturbation theory

`bhpwave` generates gravitational waveforms of extreme-mass-ratio inspirals (EMRIs) using the adiabatic approximation of black hole perturbation theory. The model is restricted to binaries in which the small body is undergoing a quasi-circular, equatorial inspiral into a rotating massive black hole. Details and applications of
the `bhpwave` model can be found on [Physical Review D](https://doi.org/10.1103/PhysRevD.109.044020), [arXiv](https://doi.org/10.48550/arXiv.2310.19706), and the Github
repository [bhpwave-article](https://github.com/znasipak/bhpwave-article). Detailed documentation is provided [here](https://bhpwave.readthedocs.io/en/latest/). 

# Installation

`bhpwave` relies on a few dependencies to install and run, namely
a C/C++ compiler (e.g., `g++`), `gsl`, `Cython`, 
`numpy`, and `python >= 3.7`, though we recommend using at least Python 3.9.
To reduce package conflicts and ensure that the proper dependencies are installed,
we recommend using Anaconda and its virtual environments.

Create a conda environment `bhpwave-env` (or whatever name you would like)
with the necessary dependencies to install `bhpwave`. For MACOSX Intel run:
```
conda create -n bhpwave-env -c conda-forge gsl Cython numpy clang_osx-64 clangxx_osx-64 python=3.9
conda activate bhpwave-env
```
This may also work for MACOSX silicon, though alternatively one should use:
For MACOSX Intel run:
```
conda create -n bhpwave-env -c conda-forge gsl Cython numpy clang_osx-arm64 clangxx_osx-arm64 python=3.9
conda activate bhpwave-env
```
See Troubleshooting.
To instead include the necessary compiler on linux run:
```
conda create -n bhpwave-env -c conda-forge gsl Cython numpy gcc_linux-64 gxx_linux-64 python=3.9
conda activate bhpwave-env
```
Next clone the :code:`bhpwave` repository from GitHub:
```
git clone https://github.com/znasipak/bhpwave.git
cd bhpwave
```
Finally, we recommend installing the package via `pip`:
```
pip install .
```

# Conda Environments with Jupyter

To run the code in a jupyter notebook, we recommend `pip` installing
the following dependencies: 
```
pip install ipykernel matplotlib
```
One can then make the environment accessible within jupyter by
running
```
python -m ipykernel install --user --name=bhpwave-env
```

# Troubleshooting

One may run into compiler issues, particularly for MACOSX silicon machines. If you are trying to resolve compiler errors, make sure to run
```
pip uninstall bhpwave
```
and to remove the `build/` and `bhpwave.egg-info` directories before trying to recompile the code. One might also need to remove `-march=native` 

# License

This package is provided under the MIT license but links to GPL-licensed libraries. 

## Licensing Consequences
> [!WARNING]
> The following is not legal advice.

To the best of the authors' knowledge, one can redistribute the code under an MIT-compatible but GPL-incompatible license only by removing the GPL-dependent libraries. Otherwise, the the entire software package (code + dependencies) falls under the GPL license.
