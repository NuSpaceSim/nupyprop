.. nuPyProp documentation master file, created by
   sphinx-quickstart on Fri Jul 16 10:47:33 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |br| raw:: html

Welcome to nuPyProp's documentation!
====================================


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`


Welcome to nuPyProp
===================

.. comment: split here

nuPyProp is a Monte Carlo simulation package which models the neutrino flux
attenuation and the energy distribution of the leptons it creates while
propagating through the Earth. The package is also a part of nuSpaceSim.




..  for computation of non-thermal radiation from
.. relativistic particle populations. It includes tools to perform MCMC fitting of
.. radiative models to X-ray, GeV, and TeV spectra using `emcee
.. <http://dan.iel.fm/emcee>`_, an affine-invariant ensemble sampler for Markov
.. Chain Monte Carlo. Naima is an `Astropy`_ affiliated
.. package.

.. comment: - Code: http://www.github.com/zblz/naima

.. comment: - Documentation: http://naima.readthedocs.org

.. comment: split here

The code contains 7 main modules:-<br/>
* :ref:`data`: This is responsible for reading and writing data to the lookup tables
hdf file and the resultant output hdf file.<br/>
* :ref:`geometry`: Mostly deals with the Earth geometry and the functions needed for
generating water and column density trajectories for propagating through the Earth.<br/>
* :ref:`cross-section`: Serves as a customizable template for the different neutrino &
anti-neutrino cross-section models.<br/>
* :ref:`energy-loss`: Also serves as a customizable template for the different charged
lepton (tau and mu) electromagnetic energy loss models.<br/>
* :ref:`propagate`: The only fortran module in the packge, which simulates the 
propagation of the particles inside the Earth. Contains some Earth geometry functions,
interpolation and propagation subroutines. Wrapped in python using f2py.<br/>
* :ref:`main`: Used for initializing parameters, cross-sections and energy loss
models. It is a pythonic interface for using the propagation module.<br/>
* :ref:`run`: Used for parsing runtime arguments and parameters for running the
code.<br/>


.. There are two main components of the package: a set of nonthermal
.. :ref:`radiative`, and a set of utility functions that make it easier to fit a
.. given model to observed spectral data (see :ref:`MCMC`).

.. Nonthermal radiative models are available for Synchrotron, inverse Compton,
.. Bremsstrahlung, and neutral pion decay processes. All of the models allow the
.. use of an arbitrary shape of the particle energy distribution, and several
.. functional models are also available to be used as particle distribution
.. functions. See :ref:`radiative` for a detailed explanation of these.

Use the sidebar on the left to access the documentation.


.. toctree::
   :hidden:
   :caption: Contents:
   :maxdepth: 2
 
   installation.rst
   data.rst
   geometry.rst
   cross-section.rst
   energy-loss.rst
   propagate.rst
   main.rst
   run.rst


.. toctree::
   :hidden:
   :maxdepth: 2

   tutorial.rst
   examples.rst

License & Attribution
---------------------

Naima is released under a 3-clause BSD style license - see the
`LICENSE.rst <https://github.com/zblz/naima/blob/master/LICENSE.rst>`_ for
details.

If you find Naima useful in your research, you can cite `Zabalza (2015)
<http://arxiv.org/abs/1509.03319>`_ (`arXiv <http://arxiv.org/abs/1509.03319>`_,
`ADS <http://adsabs.harvard.edu/abs/2015arXiv150903319Z>`_) to acknowledge its
use. The BibTeX entry for the paper is::

    @ARTICLE{naima,
       author = {{Zabalza}, V.},
        title = {naima: a Python package for inference of relativistic particle
                 energy distributions from observed nonthermal spectra},
         year = 2015,
      journal = {Proc.~of International Cosmic Ray Conference 2015},
        pages = "922",
       eprint = {1509.03319},
       adsurl = {http://adsabs.harvard.edu/abs/2015arXiv150903319Z},
    }


Contributing
------------

Please report any issues with the package `here
<https://github.com/zblz/naima/issues>`_.

All development of Naima is done through the `github repository`_, and
contributions to the code are welcome.  The development model is similar to that
of `astropy`_, so you can check the `astropy Developer Documentation
<https://astropy.readthedocs.org/en/latest/#developer-documentation>`_ if you
need information on how to make a code contribution.

.. _github repository: https://github.com/zblz/naima