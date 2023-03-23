nuPyProp Tutorial
=================

A tutorial for understanding the input & output data for the `nuPyProp <https://github.com/NuSpaceSim/nupyprop>`_ code.

While the nuPyProp `documentation <https://nupyprop.readthedocs.io/en/latest/>`_ is still work in progress,
this repository will help you visualize and understand how to plot the output data from the output file generated
through the code and will also walk you through to create user-defined models for neutrino/anti-neutrino cross-sections
charged lepton photonuclear energy loss models.

There are three jupyter notebooks in this tutorial, where ``charge_lepton_energy_loss.ipynb`` and ``neutrino_cross_sections.ipynb`` notebooks shows the user how to visualize different electromagnetic energy loss models and neutrino cross section models used in the nuPyProp code. Third notebook, ``output_plots.ipynb``, shows how to visualize the output h5 files. The output files are in Hierarchical Data Format (HDF). One way to open the h5 files is by using ViTables (https://vitables.org/Download/).

Dependencies
------------

Other than the dependencies already installed by nuPyProp, you will require:

1. jupyter-lab - 
   ``pip install jupyterlab`` or ``conda install -c conda-forge jupyterlab``
2. seaborn - 
   ``pip install seaborn`` or ``conda install seaborn``

Note: It is recommended to run the ipynb notebooks in this tutorial repository in the nuPyProp environment.

Installation
------------

1. Install nuPyProp through one of the many options outlined in the `nuPyProp repository <https://github.com/NuSpaceSim/nupyprop>`_.
2. Install the external dependencies mentioned above.
3. git clone the source code:
   ``git clone https://github.com/sameer-patel-1/nupyprop_tutorial.git``
4. Either place your nuPyProp output hdf5 file in ``plotting_tutorial/data`` or use the sample data file provided with this repo.
5. ``cd nupyprop_tutorial``
6. ``jupyter-lab *.ipynb``, where * is any of the 3 ipython notebooks in the repo.

Notes
-----

1. The ipynb notebooks in this tutorial repository should be run in the nuPyProp environment.
2. There are two output files using allm electromagnetic energy loss method in ``plotting_tutorial/data``, named ``output_nu_tau_4km_ct18nlo_allm_stochastic_1e7.h5`` and, ``output_nu_tau_4km_ct18nlo_allm_stochastic_1e8.h5`` . And another output file using bdhm electromagnetic energy loss method in ``plotting_tutorial/data``, named ``output_nu_tau_4km_ct18nlo_bdhm_stochastic_1e8.h5`` .
3. The file ``output_nu_tau_4km_ct18nlo_allm_stochastic_1e8.h5`` is used in the ``output_plots.ipynb`` notebook to show how to plot the output data from nuPyProp. This file has tau-neutrino energy ranging from 1e6 to 1e12 GeV in quarter decades. This output file has two tables, Exit_Probability and Clep_out_cdf. We didn't include the Clep_out_energy because that made the h5 file really big in size. 
4. Therefore, section 2 can't be implemented using the given example h5 files. But an example code is given to the user to visualize the reuslts, if they have Clep_out_energy table in the h5 file (which can be included by running the nuPyProp simulation code with --el flag as yes, look at nuPyProp README document for more information).
5. The file ``output_nu_tau_4km_ct18nlo_allm_stochastic_1e7.h5`` and ``output_nu_tau_4km_ct18nlo_bdhm_stochastic_1e8.h5`` has three tables, Exit_Probability, Avg_Polarization and Clep_out_cdf. This is an example file for plotting average polarization of the exiting lepton. This is shown in section 6 of the ``output_plots.ipynb`` notebook. 

Creating Custom Models
----------------------

The custom/user-defined models will have to be of the same formats as shown in the iPython notebooks
and the user should utilize the python code (`models.py <https://github.com/NuSpaceSim/nupyprop/blob/main/src/nupyprop/models/models.py>`_) which provides examples and templates on how to generate your own models. This is further more explained in the paper: "Neutrino propagation in the Earth and emerging charged leptons with nuPyProp", D. Garg, S.Patel et al. (NuSpaceSim Collaboration), `e-Print: arXiv:2209.15581 [astro-ph.HE, hep-ph] <https://doi.org/10.48550/arXiv.2209.15581>`__ 


.. |alpha|  unicode:: U+003B1 .. GREEK SMALL LETTER ALPHA
.. |beta|   unicode:: U+003B2 .. GREEK SMALL LETTER BETA
.. |chi|    unicode:: U+003C7 .. GREEK SMALL LETTER CHI
.. |Delta|  unicode:: U+00394 .. GREEK CAPITAL LETTER DELTA
.. |delta|  unicode:: U+003B4 .. GREEK SMALL LETTER DELTA
.. |epsi|   unicode:: U+003F5 .. GREEK LUNATE EPSILON SYMBOL
.. |epsis|  unicode:: U+003F5 .. GREEK LUNATE EPSILON SYMBOL
.. |epsiv|  unicode:: U+003B5 .. GREEK SMALL LETTER EPSILON
.. |eta|    unicode:: U+003B7 .. GREEK SMALL LETTER ETA
.. |Gamma|  unicode:: U+00393 .. GREEK CAPITAL LETTER GAMMA
.. |gamma|  unicode:: U+003B3 .. GREEK SMALL LETTER GAMMA
.. |Gammad| unicode:: U+003DC .. GREEK LETTER DIGAMMA
.. |gammad| unicode:: U+003DD .. GREEK SMALL LETTER DIGAMMA
.. |iota|   unicode:: U+003B9 .. GREEK SMALL LETTER IOTA
.. |kappa|  unicode:: U+003BA .. GREEK SMALL LETTER KAPPA
.. |kappav| unicode:: U+003F0 .. GREEK KAPPA SYMBOL
.. |Lambda| unicode:: U+0039B .. GREEK CAPITAL LETTER LAMDA
.. |lambda| unicode:: U+003BB .. GREEK SMALL LETTER LAMDA
.. |mu|     unicode:: U+003BC .. GREEK SMALL LETTER MU
.. |nu|     unicode:: U+003BD .. GREEK SMALL LETTER NU
.. |Omega|  unicode:: U+003A9 .. GREEK CAPITAL LETTER OMEGA
.. |omega|  unicode:: U+003C9 .. GREEK SMALL LETTER OMEGA
.. |Phi|    unicode:: U+003A6 .. GREEK CAPITAL LETTER PHI
.. |phi|    unicode:: U+003D5 .. GREEK PHI SYMBOL
.. |phis|   unicode:: U+003D5 .. GREEK PHI SYMBOL
.. |phiv|   unicode:: U+003C6 .. GREEK SMALL LETTER PHI
.. |Pi|     unicode:: U+003A0 .. GREEK CAPITAL LETTER PI
.. |pi|     unicode:: U+003C0 .. GREEK SMALL LETTER PI
.. |piv|    unicode:: U+003D6 .. GREEK PI SYMBOL
.. |Psi|    unicode:: U+003A8 .. GREEK CAPITAL LETTER PSI
.. |psi|    unicode:: U+003C8 .. GREEK SMALL LETTER PSI
.. |rho|    unicode:: U+003C1 .. GREEK SMALL LETTER RHO
.. |rhov|   unicode:: U+003F1 .. GREEK RHO SYMBOL
.. |Sigma|  unicode:: U+003A3 .. GREEK CAPITAL LETTER SIGMA
.. |sigma|  unicode:: U+003C3 .. GREEK SMALL LETTER SIGMA
.. |sigmav| unicode:: U+003C2 .. GREEK SMALL LETTER FINAL SIGMA
.. |tau|    unicode:: U+003C4 .. GREEK SMALL LETTER TAU
.. |Theta|  unicode:: U+00398 .. GREEK CAPITAL LETTER THETA
.. |theta|  unicode:: U+003B8 .. GREEK SMALL LETTER THETA
.. |thetas| unicode:: U+003B8 .. GREEK SMALL LETTER THETA
.. |thetav| unicode:: U+003D1 .. GREEK THETA SYMBOL
.. |Upsi|   unicode:: U+003D2 .. GREEK UPSILON WITH HOOK SYMBOL
.. |upsi|   unicode:: U+003C5 .. GREEK SMALL LETTER UPSILON
.. |Xi|     unicode:: U+0039E .. GREEK CAPITAL LETTER XI
.. |xi|     unicode:: U+003BE .. GREEK SMALL LETTER XI
.. |zeta|   unicode:: U+003B6 .. GREEK SMALL LETTER ZETA

