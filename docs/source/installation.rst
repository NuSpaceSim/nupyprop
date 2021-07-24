.. _installation:

Installation
============

The Nupyprop package relies heavily on compiled fortran, parallelized with openMP to
propagate neutrinos. Users however are not required to install these arcane components
in order to use nupyprop. Instead the nuSpaceSim team provides pre-compiled
implementations of these software libraries, along with their necessary dependencies.
These pre-compiled packages are available via pip and conda.

Nupyprop is compatable with versions of ``python >= 3.7``. Pip or Anaconda are needed to
install these pre-compiled versions on the user's system.

Requirements
------------

Nupyprop is a python package compiled for MacOS and Linux operating systems running on
``x86_64`` processors. User systems are required to have a working python3 interpreter
running ``python 3.7`` or later, as well as a python package manager application such as
`pip`_ or `conda`_. All other runtime requirements will be installed by the user's package
manager.

.. _pip: https://pip.pypa.io/en/stable/
.. _conda: https://docs.conda.io/en/latest/miniconda.html

Install with pip
----------------

Installing nupyprop into the active python3 instance is easily accomplished with:

::

  python3 -m pip install nupyprop

If the user does not have access to this space, because of system permissions issues or
other related errors, nupyprop can be installed locally into the user's home directory
via:

::

  python3 -m pip install --user nupyprop

Check out the `nupyprop PyPi page <https://pypi.org/project/nupyprop/>`_.


Install with conda
------------------

We recommend installing nupyprop into a conda environment like so. In
this example the name of the environment is ``nuspacesim``

::

   conda create -n nuspacesim -c conda-forge -c nuspacesim nupyprop
   conda activate nuspacesim

To install the nupyprop package into an existing conda environment named ``nuspacesim``,
execute the following:

::

  conda install -n nuspacesim -c conda-forge -c nuspacesim nupyprop


Install from source
-------------------

Nupyprop is compliant with both :PEP:`517` and :PEP:`518`, and relies on pip for as a
PEP 517 compatable build tool to build the package into a binary wheel and install it on
host systems.

Requirements
~~~~~~~~~~~~

If you wish to install the code from source, you will need a working fortran compiler such
as gfortran, and an OpenMP implementation. On MacOS and Ubuntu both are provided as part
of the gcc compiler suite.

You will also need ``pip`` as all additional python requirements such as ``numpy`` will
be installed by pip.

.. tip::
    Install gcc on your system using a package manager.

    On Ubuntu Linux with apt-get use:

    ::

       sudo apt-get install gcc gfortran openmp

    On MacOS with `Homebrew`_ use:

    ::

      brew install gcc

.. _Homebrew: https://brew.sh/

To install from source:

1. Obtain the source code from the nupyprop git repository:
   ::

      git clone https://github.com/NuSpaceSim/nupyprop
      cd nupyprop/
2. Use pip install to build and install the package:
   ::

      python3 -m pip install .



