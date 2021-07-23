nupyprop
========

Propagate neutrinos through the earth.

A python package and command line utility, including fortran for
performance with openMP.

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

**Example** for running tau propagation for 107 GeV neutrinos at 10
degrees with a statistics of 107 particles with stochastic energy loss &
with all other parameters as defaults:

``nupyprop -e 7 -a 10 -t stochastic -s 1e7``

**Run parameters** are defined in run.py. Different switches are
described as follows:

1. -e or –energy: incoming neutrino energy in log_10(GeV). Works for
   single energy or multiple energies. For multiple energies, separate
   energies with commas eg. 7,8,9,10,11. Default energies are
   107,107.25,107.5,…1011 GeV.

2. -a or –angle: slant Earth angles in degrees. Works for single angle
   or multiple angles. For multiple angles, separate angles with commas
   eg. 1,3,5,7,10. Default angles are 1->35 degrees, in steps of 1
   degree.

3. -i or –idepth: depth of ice/water in km. Default value is 4 km.

4. -l or –lepton: flavor of lepton used to propagate. Can be either muon
   or tau. Default is tau.

5. -n or –nu_type: type of neutrino matter. Can be either neutrino or
   anti-neutrino. Default is Neutrino.

6. -t or –energy_loss: energy loss type for lepton - can be stochastic
   or continuous. Default is stochastic.

..    ~~7. -m or --material: material used in electromagnetic energy loss; not used in main program, only used for running energy_loss.py individually. Default is rock.~~

7.  -x or –xc_model: neutrino/anti-neutrino cross-section model used.
    For a list of model names, see
    lookup_tables.h5/Neutrino_Cross_Sections. Default is ct18nlo.

8.  -p or –pn_model: photonuclear interaction energy loss model used.
    For now, can be either bb (Bezrukov-Bugaev), allm (Abramowicz,
    Levin, Levy, Maor), bdhm (Block, Durand, Ha, McKay) or ckmt
    (Capella, Kaidalov, Merino, Tran). Default is allm.

9.  -f or –fac_nu: rescaling factor for SM cross-sections. Default is 1.

10. -s or –stat: statistics. Default is 1e7 particles.

[STRIKEOUT:WARNING: Running the code will replace the results of the
output.h5 file. Backing up previous output files is recommended (or just
ask me for a pre-populated output file if need be, for now, since the
pre-populated output file is >25 MB). Future fixes include overwrite
warnings for the user.] Fixed!

**Viewing output results**: output_*.h5 will contain the results of the
code after it has finished running. In the terminal, run ``vitables``
and open the output_*.h5 file to view the output results.

output_*.h5 naming convention is as follows: **output_A_B_C_D_E_F_G**,
where

A = Neutrino type: nu is for neutrino & anu is for anti-neutrino. B =
Lepton_type: tau is for tau leptons & muon is for muons. C = idepth:
depth of water layer in km. D = Neutrino (or anti-neutrino)
cross-section model. E = Lepton photonuclear energy loss model. F =
Energy loss type: can be stochastic or continuous. G = Statistics (ie.
no. of neutrinos/anti-neutrinos injected).

.. table:: **Code Execution Timing Table for Taus**:

   ============== ================ ====================== ====== =================== ==========
   Charged Lepton Energy Loss Type E<sub>&nu;</sub> [GeV] Angles N<sub>&nu;;in</sub> Time (hrs)
   ============== ================ ====================== ====== =================== ==========
   :math:`{tau}`  Stochastic       10\ :sup:`7`           1-35   10\ :sup:`8`        1.07*, 0.26***  
   :math:`{tau}`  Continuous       10\ :sup:`7`           1-35   10\ :sup:`8`        0.88*           
   :math:`{tau}`  Stochastic       10\ :sup:`8`           1-35   10\ :sup:`8`        6.18*, 1.53***  
   :math:`{tau}`  Continuous       10\ :sup:`8`           1-35   10\ :sup:`8`        5.51*           
   :math:`{tau}`  Stochastic       10\ :sup:`9`           1-35   10\ :sup:`8`        27.96*, 5.08*** 
   :math:`{tau}`  Continuous       10\ :sup:`9`           1-35   10\ :sup:`8`        19.11*          
   :math:`{tau}`  Stochastic       10\ :sup:`10`          1-35   10\ :sup:`8`        49.80*, 12.43***
   :math:`{tau}`  Continuous       10\ :sup:`10`          1-35   10\ :sup:`8`        35.59*          
   :math:`{tau}`  Stochastic       10\ :sup:`11`          1-35   10\ :sup:`8`        12.73***        
   :math:`{tau}`  Continuous       10\ :sup:`11`          1-35   10\ :sup:`8`        -               
   ============== ================ ====================== ====== =================== ==========


.. table:: **Code Execution Timing Table for Muons**:

  ============== ================ ====================== ================================= =================== ==========
  Charged Lepton Energy Loss Type E<sub>&nu;</sub> [GeV] Angles                            N<sub>&nu;;in</sub> Time (hrs)
  ============== ================ ====================== ================================= =================== ==========
  :math:`{mu}`   Stochastic       10\ :sup:`6`           1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:8`         -        
  :math:`{mu}`   Continuous       10\ :sup:`6`           1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:8`         0.95*    
  :math:`{mu}`   Stochastic       10\ :sup:`7`           1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:8`         -        
  :math:`{mu}`   Continuous       10\ :sup:`7`           1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:8`         3.19*    
  :math:`{mu}`   Stochastic       10\ :sup:`8`           1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:8`         -        
  :math:`{mu}`   Continuous       10\ :sup:`8`           1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:8`         5.17*    
  :math:`{mu}`   Stochastic       10\ :sup:`9`           1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:8`         111.77** 
  :math:`{mu}`   Continuous       10\ :sup:`9`           1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:8`         7.42*    
  :math:`{mu}`   Stochastic       10\ :sup:`10`          1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:8`         98.17*   
  :math:`{mu}`   Continuous       10\ :sup:`10`          1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:8`         9.76*    
  :math:`{mu}`   Stochastic       10\ :sup:`11`          1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:8`         -        
  :math:`{mu}`   Continuous       10\ :sup:`11`          1,2,3,5,7,10,12,15,17,20,25,30,35 10\ :sup:8`         -        
  ============== ================ ====================== ================================= =================== ==========

\* - Intel Core i7-8750H; 6 cores & 12 threads. \*\* - Intel Core
i5-10210; 4 cores & 8 threads. \**\* - UIowa Argon cluster; 56 cores.

**For debugging/development:** The correct order to look at the code is
in the following order:

1. *data.py*: contains the reading/writing modules from/to the hdf5
   files.
2. *geometry_py.py*: contains the Earth geometry modules (including
   PREM) for use with python.
3. *cross_section.py*: contains neutrino/anti-neutrino cross_section
   models.
4. *energy_loss.py*: contains lepton energy loss models.
5. *propagate.f90*: heart of the code; contains fortran modules to
   interpolate between geometry variables, cross-sections, energy loss
   parameters & propagate neutrinos and leptons through the Earth.
6. *main.py*: forms the main skeleton of the code; propagates the
   neutrinos and leptons, and calculates the p_exit and collects
   outgoing lepton energies.
7. *run.py*: contains all the run parameters and variables needed for
   all the other .py files.

.. figure:: /figures/nupyprop_uml_full.png
   :alt: UML Diagram

   UML Diagram

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
