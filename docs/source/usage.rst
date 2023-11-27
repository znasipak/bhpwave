Getting Started
===============

**BHPWave** is designed to be a simple waveform generator for
small bodies undergoing quasi-circular inspirals into rotating massive
black holes (MBHs).

.. note::
   Installation has only been tested on MACOSX (Intel chip) and linux
   machines. We are working on updating installation instructions for
   MACOSX (M1 chip) and Windows users.

Installation
------------

BHPWave relies on a few dependencies to install and run, namely
a C/C++ compiler (e.g., :code:`g++`), :code:`gsl`, :code:`Cython`, 
:code:`numpy`, and :code:`python >= 3.7`.
To reduce package conflicts and ensure that the proper dependencies are installed,
we recommend using Anaconda and its virtual environments.

Create a conda environment :code:`bhpwave-env` (or whatever name you would like)
with the necessary dependencies to install :code:`bhpwave`:

.. code-block:: console

    conda create -n bhpwave-env -c conda-forge gsl Cython numpy python=3.7
    conda activate bhpwave-env

To include the necessary compiler on MACOSX run:

.. code-block:: console

    conda install clang_osx-64 clangxx_osx-64

To include the necessary compiler on linux run:

.. code-block:: console

    conda install gcc_linux-64 gxx_linux-64

Next clone the :code:`bhpwave` repository from GitHub:

.. code-block:: console

    git clone https://github.com/znasipak/bhpwave.git
    cd bhpwave

Finally, we recommend installing the package via :code:`pip`:

.. code-block:: console

    pip install .

Conda Environments with Jupyter
-------------------------------

To run the code in a jupyter notebook, we recommend :code:`pip` installing
the following dependencies: 

.. code-block:: console

    pip install ipykernel matplotlib

One can then make the environment accessible within jupyter by
running

.. code-block:: console

    python -m ipykernel install --user --name=bhpwave-env

Uninstalling
------------

If the package is installed using :code:`pip`, then one can easily uninstall the package
using :code:`pip` as well

.. code-block:: console

    pip uninstall bhpwave

To clean the repository, one will also need to remove the directories
:code:`build` and :code:`bhpwave.egg-info` from the main repository

