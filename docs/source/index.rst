.. nuPyProp documentation master file, created by
   sphinx-quickstart on Fri Jul 16 10:47:33 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |br| raw:: html

   <br />


Welcome to nuPyProp
===================

nuPyProp is a Monte Carlo simulation package which models the neutrino flux
attenuation and the energy distribution of the leptons it creates while
propagating through the Earth. The package is also a part of nuSpaceSim.


The code package contains 7 main modules:-

1. :ref:`data`: This is responsible for reading and writing data to the lookup tables
hdf file and the resultant output hdf file.

2. :ref:`geometry`: Mostly deals with the Earth geometry and the functions needed for
generating water and column density trajectories for propagating through the Earth.

3. :ref:`cross-section`: Serves as a customizable template for the different neutrino &
anti-neutrino cross-section models.

4. :ref:`energy-loss`: Also serves as a customizable template for the different charged
lepton (tau and mu) electromagnetic energy loss models.

5. :ref:`propagate`: The only fortran module in the packge, which simulates the 
propagation of the particles inside the Earth. It contains some Earth geometry functions,
interpolation and propagation subroutines, and has been wrapped in python using f2py.

6. :ref:`main`: Used for initializing parameters, cross-sections and energy loss
models. It is a pythonic interface for using the propagation module.

7. :ref:`run`: Used for parsing runtime arguments and parameters for running the
code.

Use the sidebar on the left to access the documentation for each module.

.. toctree::
   :hidden:
   :caption: Contents:
   :maxdepth: 2
 
   installation.rst
   run.rst
   main.rst
   data.rst
   geometry.rst
   cross-section.rst
   energy-loss.rst
   propagate.rst

.. toctree::
   :hidden:
   :maxdepth: 2

   tutorial.rst
   examples.rst
   development.rst
