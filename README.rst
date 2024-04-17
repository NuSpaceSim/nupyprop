nuPyProp
========

Propagate neutrinos through the earth.

A python package and command line utility, including fortran for
performance with openMP.

Documentation (WIP): https://nupyprop.readthedocs.io/en/latest/

**Note:** While the documentation is currently WIP, users and developers should consult the 
nuPyProp tutorial folder called **tutorial** in the current repository, for visualizing output from the code and creating user-defined models.

**Citation:** Please cite "Neutrino propagation in the Earth and emerging charged leptons with nuPyProp", Diksha Garg et al JCAP01(2023)041. DOI: 10.1088/1475-7516/2023/01/041 (`Paper link <https://iopscience.iop.org/article/10.1088/1475-7516/2023/01/041>`__), published in Journal of Cosmology and Astroparticle Physics.

**Acknowledgments:** This work is supported by NASA grants 80NSSC19K0626 at the University of Maryland,Baltimore County, 80NSSC19K0460 at the Colorado School of Mines, 80NSSC19K0484 at theUniversity of Iowa, and 80NSSC19K0485 at the University of Utah, 80NSSC18K0464 at LehmanCollege, and under proposal 17-APRA17-0066 at NASA/GSFC and JPL.


Installation
------------

with pip
~~~~~~~~

::

   python3 -m pip install nupyprop

with conda
~~~~~~~~~~

We recommend installing nupyprop into a conda environment like so. In
this example the name of the environment is “nupyprop”

::

   conda create -n nupyprop -c conda-forge -c nuspacesim nupyprop
   conda activate nupyprop

Usage
-----

``nupyprop --help``

**Example** for running tau propagation for 10\ :sup:`7` GeV neutrinos at 10
degrees with 10\ :sup:`7` neutrinos injected with stochastic energy loss &
with all other parameters as defaults:

``nupyprop -e 7 -a 10 -t stochastic -s 1e7``

**Run parameters** are defined in run.py. Different switches are
described as follows:

1. ``-e`` or ``--energy``: incoming neutrino energy in log_10(GeV). Works for
   single energy or multiple energies. For multiple energies, separate
   energies with commas eg. 6,7,8,9,10,11. Default energies are
   10\ :sup:`6`,10\ :sup:`6.25`,10\ :sup:`6.5`,…10\ :sup:`11` GeV.

2. ``-a`` or ``--angle``: slant Earth angles in degrees. Works for single angle
   or multiple angles. For multiple angles, separate angles with commas
   eg. 1,3,5,7,10. Default angles are 1->42 degrees, in steps of 1
   degree.

3. ``-i`` or ``--idepth``: depth of ice/water in km. Default value is 4 km.

4. ``-cl`` or ``--charged_lepton``: flavor of charged lepton used to propagate. Can be either
   muon or tau. Default is tau.

5. ``-n`` or ``--nu_type``: type of neutrino matter. Can be either neutrino or
   anti-neutrino. Default is neutrino.

6. ``-t`` or ``--energy_loss``: energy loss type for lepton - can be stochastic
   or continuous. Default is stochastic.

7.  ``-x`` or ``--xc_model``: neutrino/anti-neutrino cross-section model used.
    Can be from the pre-defined set of models (see xc-table_) or custom.
    Default is ct18nlo.

8.  ``-p`` or ``--pn_model``: photonuclear interaction energy loss model used.
    Can be from the pre-defined set of models (see pn-table_) or custom.
    Default is allm.

9.  ``-el`` or ``--energy_lepton``: option to print exiting charged lepton's final energy in
    output file.
    Default is no

10.  ``-f`` or ``--fac_nu``: rescaling factor for BSM cross-sections. Default is 1.

11. ``-s`` or ``--stats``: statistics or no. of injected neutrinos. Default is 1e7
    neutrinos.
    
12. ``-htc`` or ``--htc_mode``: High throughput computing mode. If set to yes,
    the code will be optimized to run in high throughput computing mode.
    Default is no.
    
**Note**: This program uses OpenMP for propagating the huge number of neutrinos injected.
For this purpose, the code will use **all** the threads available in your processor by default.
To control the number of threads used for running the code, use ``export OMP_NUM_THREADS=x``, 
where ``x`` is the number of threads you want the code to run with.

**Viewing output results**: output_*.h5 will contain the results of the
code after it has finished running. In the terminal, run ``vitables``
(optional dependency) and open the output_*.h5 file to view the output
results.

output_*.h5 naming convention is as follows: **output_A_B_Ckm_D_E_F_G**,
where

| A = Neutrino type: nu is for neutrino & anu is for anti-neutrino.
| B = Charged lepton type: tau is for tau leptons & muon is for muons.
| C = idepth: depth of water layer (in km).
| D = Neutrino (or anti-neutrino) cross-section model.
| E = Charged lepton photonuclear energy loss model.
| F = Energy loss type: can be stochastic or continuous.
| G = Statistics (ie. no. of neutrinos/anti-neutrinos injected).

Model Tables
------------

.. _xc-table:
   
   +--------------------------------------------+--------------------------------------------------------------------------------------------------+
   | Neutrino/Anti-Neutrino Cross-Section Model |                                             Reference                                            |
   +============================================+==================================================================================================+
   |    Abramowicz, Levin, Levy, Maor (ALLM)    | `hep-ph/9712415 <https://arxiv.org/abs/hep-ph/9712415>`_,                                        |
   |                                            | `Phys. Rev. D 96, 043003 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.96.043003>`_    |
   +--------------------------------------------+--------------------------------------------------------------------------------------------------+
   |       Block, Durand, Ha, McKay (BDHM)      | `Phys. Rev. D 89, 094027 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.89.094027>`_,   |
   |                                            | `Phys. Rev. D 96, 043003 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.96.043003>`_    |
   +--------------------------------------------+--------------------------------------------------------------------------------------------------+
   |                 CTEQ18-NLO                 | `Phys. Rev. D 103, 014013 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.103.014013>`_, |
   |                                            | `Phys. Rev. D 81, 114012 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.81.114012>`_    |
   +--------------------------------------------+--------------------------------------------------------------------------------------------------+
   |       Connolly, Thorne, Waters (CTW)       | `Phys. Rev. D 83, 113009 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.83.113009>`_    |
   +--------------------------------------------+--------------------------------------------------------------------------------------------------+
   |                   nCTEQ15                  | `Phys. Rev. D 93, 085037 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.93.085037>`_,   |
   |                                            | `Phys. Rev. D 81, 114012 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.81.114012>`_    |
   +--------------------------------------------+--------------------------------------------------------------------------------------------------+
   |                User Defined                | See `nuPyProp tutorial repository <https://research-git.uiowa.edu/spatel31/nupyprop_tutorial>`__ |
   +--------------------------------------------+--------------------------------------------------------------------------------------------------+
   


.. _pn-table:

   +-----------------------------------------------+-------------------------------------------------------------------------------------------------+
   | Charged Lepton Photonuclear Energy Loss Model |                                            Reference                                            |
   +===============================================+=================================================================================================+
   |      Abramowicz, Levin, Levy, Maor (ALLM)     | `hep-ph/9712415 <https://arxiv.org/abs/hep-ph/9712415>`_,                                       |
   |                                               | `Phys. Rev. D 63, 094020 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.63.094020>`_   |
   +-----------------------------------------------+-------------------------------------------------------------------------------------------------+
   |              Bezrukov-Bugaev (BB)             | `Yad. Fiz. 33, 1195 <https://inspirehep.net/literature/170124>`_,                               |
   |                                               | `Phys. Rev. D 63, 094020 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.63.094020>`_   |
   +-----------------------------------------------+-------------------------------------------------------------------------------------------------+
   |        Block, Durand, Ha, McKay (BDHM)        | `Phys. Rev. D 89, 094027 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.89.094027>`_,  |
   |                                               | `Phys. Rev. D 63, 094020 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.63.094020>`_   |
   +-----------------------------------------------+-------------------------------------------------------------------------------------------------+
   |     Capella, Kaidalov, Merino, Tran (CKMT)    | `Eur. Phys. J. C 10, 153 <https://arxiv.org/abs/hep-ph/9806367>`_                               |
   |                                               | `Phys. Rev. D 63, 094020 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.63.094020>`_   |
   +-----------------------------------------------+-------------------------------------------------------------------------------------------------+
   |                  User Defined                 | See `nuPyProp tutorial repository <https://research-git.uiowa.edu/spatel31/nupyprop_tutorial>`__|
   +-----------------------------------------------+-------------------------------------------------------------------------------------------------+

Code Execution Timing Tables
----------------------------
.. _tau-table:

   ============== ================ ==================== ====== =================== ===============
   Charged Lepton Energy Loss Type E\ :sub:`|nu|` [GeV] Angles N\ :sub:`|nu|;;in`    Time (hrs)
   ============== ================ ==================== ====== =================== ===============
   |tau|          Stochastic       10\ :sup:`7`         1-35   10\ :sup:`8`        1.07*, 0.26***  
   |tau|          Continuous       10\ :sup:`7`         1-35   10\ :sup:`8`        0.88*           
   |tau|          Stochastic       10\ :sup:`8`         1-35   10\ :sup:`8`        6.18*, 1.53***  
   |tau|          Continuous       10\ :sup:`8`         1-35   10\ :sup:`8`        5.51*           
   |tau|          Stochastic       10\ :sup:`9`         1-35   10\ :sup:`8`        27.96*, 5.08*** 
   |tau|          Continuous       10\ :sup:`9`         1-35   10\ :sup:`8`        19.11*          
   |tau|          Stochastic       10\ :sup:`10`        1-35   10\ :sup:`8`        49.80*, 12.43***
   |tau|          Continuous       10\ :sup:`10`        1-35   10\ :sup:`8`        35.59*          
   |tau|          Stochastic       10\ :sup:`11`        1-35   10\ :sup:`8`        12.73***        
   |tau|          Continuous       10\ :sup:`11`        1-35   10\ :sup:`8`        -               
   ============== ================ ==================== ====== =================== ===============


.. _mu-table:

  ============== ================ ==================== ================================= ================== ==========
  Charged Lepton Energy Loss Type E\ :sub:`|nu|` [GeV] Angles                            N\ :sub:`|nu|;;in` Time (hrs)
  ============== ================ ==================== ================================= ================== ==========
  |mu|           Stochastic       10\ :sup:`6`         1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:`8`       -        
  |mu|           Continuous       10\ :sup:`6`         1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:`8`       0.95*    
  |mu|           Stochastic       10\ :sup:`7`         1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:`8`       -        
  |mu|           Continuous       10\ :sup:`7`         1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:`8`       3.19*    
  |mu|           Stochastic       10\ :sup:`8`         1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:`8`       -        
  |mu|           Continuous       10\ :sup:`8`         1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:`8`       5.17*    
  |mu|           Stochastic       10\ :sup:`9`         1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:`8`       111.77** 
  |mu|           Continuous       10\ :sup:`9`         1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:`8`       7.42*    
  |mu|           Stochastic       10\ :sup:`10`        1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:`8`       98.17*   
  |mu|           Continuous       10\ :sup:`10`        1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:`8`       9.76*    
  |mu|           Stochastic       10\ :sup:`11`        1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:`8`       -        
  |mu|           Continuous       10\ :sup:`11`        1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:`8`       -        
  ============== ================ ==================== ================================= ================== ==========

\* - Intel Core i7-8750H; 6 cores & 12 threads. \*\* - Intel Core
i5-10210; 4 cores & 8 threads. \**\* - UIowa Argon cluster; 56 cores.

**For debugging/development:** The correct order to look at the code is
in the following order:

1. *data.py*: contains functions for reading/writing from/to hdf5 files.
2. *geometry.py*: contains the Earth geometry modules (including
   PREM) for use with python/fortran.
3. *models.py*: contains neutrino cross-section & charged lepton energy loss model templates.
4. *propagate.f90*: heart of the code; contains fortran modules to
   interpolate between geometry variables, cross-sections, energy loss
   parameters & propagate neutrinos and charged leptons through the Earth.
5. *main.py*: forms the main skeleton of the code; propagates the
   neutrinos and charged leptons, and calculates the p_exit and collects
   outgoing lepton energies.
6. *run.py*: contains all the run parameters and variables needed for
   all the other .py files.

Developing the code on Ubuntu
-----------------------------

These notes should help developers of this code build and install the
package locally using a pep518 compliant build system (pip).

1. Install the non-pypi required dependencies as described for users
   above.
2. Install a fortran compiler. ex: ``sudo apt-get install gfortran``
3. git clone the source code:
   ``git clone git@github.com:NuSpaceSim/nupyprop.git``
4. ``cd nupyprop``
5. build and install the package in ‘editable’ mode
   ``python3 -m pip install -e .``

Developing the code on MacOS
----------------------------

These notes should help developers of this code build and install the
package locally using a pep518 compliant build system (pip). *Currently
we do not support the default system python3 on MacOS* which is out of
date and missing critical functionality. Use the homebrew python
instead, or a ``virtualenv``, or a conda environment.

1. Install the non-pypi required dependencies as described for users
   above.
2. Install a fortran compiler. ex: ``brew install gcc``
3. git clone the source code:
   ``git clone git@github.com:NuSpaceSim/nupyprop.git``
4. ``cd nupyprop``
5. build and install the package in ‘editable’ mode
   ``python3 -m pip install -e .``

.. This data file has been placed in the public domain.
.. Derived from the Unicode character mappings available from
   <http://www.w3.org/2003/entities/xml/>.
   Processed by unicode2rstsubs.py, part of Docutils:
   <http://docutils.sourceforge.net>.

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
