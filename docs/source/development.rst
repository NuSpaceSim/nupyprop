.. _development:

Development
===========

The following includes note and information for nupyprop developers to help orient
themselves with respect to the source code and the build system.

Developing with the source code
-------------------------------

The source code is available at https://github.com/NuSpaceSim/nupyprop. This is a git
repository hosted on the nuSpaceSim github organization. You can clone this repository
locally for version control management with

::

  git clone https://github.com/NuSpaceSim/nupyprop

Then simply ``cd nupyprop`` to enter the nupyprop directory and begin development.

Nupyprop uses a :PEP: 518 compliant, src-tree package layout.
Important `Configuration Files`_ are mostly kept at the top level of the repo, while
source code is in the ``src/nupyprop`` directory. Other important directories are:

 - ``docs``, with source code for generating documentation
 - ``figures`` for holding static images
 - ``recipe`` for building with non- :PEP: 517 compliant systems like conda
 - ``tests`` for unit tests.

A final important directory is the hidden ``.github`` directory
which holds the github Actions files for our CI pipeline.


Compiling the code for development
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Nupyprop is compliant with both :PEP:`517` and :PEP:`518`, and relies on pip for as a
PEP 517 compatable build tool to build the package into a binary wheel and install it on
host systems.

You will need a working Fortran90 compiler such as gfortran, and an OpenMP implementation.
On MacOS and Ubuntu both are provided as part of the latest gcc compiler suite.

You will also need ``pip`` as all additional python requirements such as ``numpy`` will
be installed by pip.

.. hint::
    Install gcc on your system using a package manager.

    On Ubuntu Linux with apt-get use: ``sudo apt-get install gcc gfortran openmp``

    On MacOS with `Homebrew`_ use: ``brew install gcc``

.. _Homebrew: https://brew.sh/


Developers will want to install the package in "editable" mode in order to make fast,
iterative changes to the python code without having to recompile. pip install to build
and install the package:
::

   python3 -m pip install -e .

All changes to the fortran code will require recompiling. Simply rerun the above command
to generate the new shared library.

``numpy.distutils`` and Fortran
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fortran module ``propagate`` is compiled with numpy distutils. Minimal specification included in
setup.py

Try to minimize the amount of code here so the build system can be as declarative as
reasonable.

Version management
~~~~~~~~~~~~~~~~~~

nupyprop uses setuptools_scm for versioning. This means the version number is generated
at build time from the latest tag in the git repository. To set a new version number
use git tag to tag your commit with a semantic version scheme.

::

   git tag v0.0.1

If you have added new commits but are not ready to tag a new version, setuptools_scm will
update the version by adding a ``-postN`` to the version string, where ``N`` is the number
of commits past the latest tag version your git repository is.

the version is written to a file named ``_version.py``, which is ignored by git.


Configuration Files
~~~~~~~~~~~~~~~~~~~

Build requirements specified in pyproject.toml

Metadata and Runtime requirements in setup.cfg

Tooling configuration in goth pyproject.toml and setup.cfg based on tool support.

Fortran extensions compiled in setup.py

Nupyprop as a console_script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once installed users use the command line executable. This is managed using setuptools'
entry_points console_script feature, in setup.cfg.

Running the tests
-----------------

Tests require pytest and are held in ``tests`` directory. To run all the tests locally,
run:
::

  pytest tests/

To test in an isolated environment for your python version use tox. If my python version
is 3.9 I run:
::

  tox -e py39

Building the wheel
------------------

::

  python3 -m pip wheel -w dist --no-deps .


Building the conda package
--------------------------
::

  conda install conda-build #if not installed
  conda build -c conda-forge recipe/

Building the documentation
--------------------------
::

  tox -e docs

Using the CI pipeline
---------------------

`https://github.com/NuSpaceSim/nupyprop/actions/workflows/pypi-build-test-publish.yml`

Run automatically on new pushes, pull-requests and release tags.

Can be run manually with the workflow dispatch "Run Workflow" option.
