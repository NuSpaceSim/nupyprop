### **Dependencies:**
python 3+ (anaconda recommended), numpy, pandas, scipy, sympy, matplotlib, fortran compiler (gfortran with openMP enabled), vitables (conda-forge; optional).

**Note:** For those who had already created a nuPyProp environment in anaconda before, you can remove it by:</br>`conda env remove --n nupyprop`

Once you have anaconda up and running, run `conda env create -f dependencies.yml`. This will create a new environment named 'nupyprop' and will install all the dependencies. It may take a while to do so. After it finishes installing, run `conda activate nupyprop` and you should be ready to run the code.

**Run parameters** are defined in run.py. Different switches are described as follows:

1. -e or --energy: incoming neutrino energy in GeV. Works for single energy or multiple energies. For multiple energies, separate energies with commas eg. 1e7,1e8,1e9. Default energies are 1e7,1e8,1e9,1e10,1e11 GeV.

2. -a or --angle: slant Earth angles in degrees. Works for single angle or multiple angles. For multiple angles, separate angles with commas eg. 1,3,5,7,10. Default angles are 1->35 degrees, in steps of 1 degree.

3. -i or --idepth: depth of ice/water in km. Default value is 4 km.

4. -l or --lepton: flavor of lepton used to propagate. Can be either muon or tau. Default is tau.

5. -t or --energy_loss: energy loss type for lepton - can be stochastic or continuous. Default is stochastic

5. -m or --material: material used in electromagnetic energy loss; not used in main program, only used for running energy_loss.py individually. Default is rock.

6. -x or --xc_model: neutrino/anti-neutrino cross-section model used. For a list of model names, see lookup_tables.h5/Neutrino_Cross_Sections. Default is CT18-NLO. 

7. -p or --pn_model: photonuclear interaction energy loss model used. For now, can be either bb (Bezrukov-Bugaev), allm (Abramowicz, Levin, Levy, Maor), bdhm (Block, Durand, Ha) or ckmt (Capella, Kaidalov, Merino, Tran). Default is allm.

8. -f or --fac_nu: rescaling factor for SM cross-sections. Default is 1.

9. -s or --stat: statistics. Default is 1e7 particles.

**Example** for running tau propagation for 10<sup>7</sup> GeV neutrinos at 10 degrees with a statistics of 10<sup>7</sup> particles with stochastic energy loss & with all other parameters as defaults:
1. `f2py -m propagate --fcompiler=gfortran --f90flags='-fopenmp' -lgomp -c propagate.f90` This is to compile the fortran file. Be sure to change `--fcompiler=gfortran` to your fortran compiler in case you use a different one. In the future, this will be done through a makefile.
2. `python run.py -e 1e7 -a 10 -t stochastic -s 1e7`

**WARNING**: Running the code will replace the results of the output.h5 file. Backing up previous output files is recommended (or just ask me for a pre-populated output file if need be, for now, since the pre-populated output file is >25 MB). Future fixes include overwrite warnings for the user.

**Viewing output results**:
output.h5 will contain the results of the code after it has finished running. In the terminal, run `vitables` and open the output.h5 file to view the output results.

**Code Execution Timing Table for Taus**:
Charged Lepton | Energy Loss Type | E<sub>&nu;</sub> [GeV] | Angles | N<sub>&nu;;in</sub> | Time (hrs) |
|---|---|---|---|---|---|
| &tau; | Stochastic | 10<sup>7</sup> | 1-35 | 10<sup>8</sup> | 1.07*, 0.26*** |
| &tau; | Continuous | 10<sup>7</sup> | 1-35 | 10<sup>8</sup> | 0.88* |
| &tau; | Stochastic | 10<sup>8</sup> | 1-35 | 10<sup>8</sup> | 6.18*, 1.53*** |
| &tau; | Continuous | 10<sup>8</sup> | 1-35 | 10<sup>8</sup> | 5.51* |
| &tau; | Stochastic | 10<sup>9</sup> | 1-35 | 10<sup>8</sup> | 27.96*, 5.08*** |
| &tau; | Continuous | 10<sup>9</sup> | 1-35 | 10<sup>8</sup> | 19.11* |
| &tau; | Stochastic | 10<sup>10</sup> | 1-35 | 10<sup>8</sup> | 49.80*, 12.43*** |
| &tau; | Continuous | 10<sup>10</sup> | 1-35 | 10<sup>8</sup> | 35.59* |
| &tau; | Stochastic | 10<sup>11</sup> | 1-35 | 10<sup>8</sup> | 12.73*** |
| &tau; | Continuous | 10<sup>11</sup> | 1-35 | 10<sup>8</sup> | - |

**Code Execution Timing Table for Muons**:
Charged Lepton | Energy Loss Type | E<sub>&nu;</sub> [GeV] | Angles | N<sub>&nu;;in</sub> | Time (hrs) |
|---|---|---|---|---|---|
| &mu; | Stochastic | 10<sup>6</sup> | 1,2,3,5,7,10,12,15,17,20,25,30,35 | 10<sup>8</sup> | - |
| &mu; | Continuous | 10<sup>6</sup> | 1,2,3,5,7,10,12,15,17,20,25,30,35 | 10<sup>8</sup> | 0.95* |
| &mu; | Stochastic | 10<sup>7</sup> | 1,2,3,5,7,10,12,15,17,20,25,30,35 | 10<sup>8</sup> | - |
| &mu; | Continuous | 10<sup>7</sup> | 1,2,3,5,7,10,12,15,17,20,25,30,35 | 10<sup>8</sup> | 3.19* |
| &mu; | Stochastic | 10<sup>8</sup> | 1,2,3,5,7,10,12,15,17,20,25,30,35 | 10<sup>8</sup> | - |
| &mu; | Continuous | 10<sup>8</sup> | 1,2,3,5,7,10,12,15,17,20,25,30,35 | 10<sup>8</sup> | 5.17* |
| &mu; | Stochastic | 10<sup>9</sup> | 1,2,3,5,7,10,12,15,17,20,25,30,35 | 10<sup>8</sup> | 111.77** |
| &mu; | Continuous | 10<sup>9</sup> | 1,2,3,5,7,10,12,15,17,20,25,30,35 | 10<sup>8</sup> | 7.42* |
| &mu; | Stochastic | 10<sup>10</sup> | 1,2,3,5,7,10,12,15,17,20,25,30,35 | 10<sup>8</sup> | 98.17* |
| &mu; | Continuous | 10<sup>10</sup> | 1,2,3,5,7,10,12,15,17,20,25,30,35 | 10<sup>8</sup> | 9.76* |
| &mu; | Stochastic | 10<sup>11</sup> | 1,2,3,5,7,10,12,15,17,20,25,30,35 | 10<sup>8</sup> | - |
| &mu; | Continuous | 10<sup>11</sup> | 1,2,3,5,7,10,12,15,17,20,25,30,35 | 10<sup>8</sup> | - |

\* - Intel Core i7-8750H; 6 cores & 12 threads.</br>
\** - Intel Core i5-10210; 4 cores & 8 threads.</br>
\*** - UIowa Argon cluster; 56 cores.

**For debugging/development:**
The correct order to look at the code is in the following order:

1. _data.py_: contains the reading/writing modules from/to the hdf5 files.
2. _geometry_py.py_: contains the Earth geometry modules (including PREM) for use with python.
3. _cross_section.py_: contains neutrino/anti-neutrino cross_section models.
4. _energy_loss.py_: contains lepton energy loss models.
5. _propagate.f90_: heart of the code; contains fortran modules to interpolate between geometry variables, cross-sections, energy loss parameters & propagate neutrinos and leptons through the Earth.
6. _main.py_: forms the main skeleton of the code; propagates the neutrinos and leptons, and calculates the p_exit and collects outgoing lepton energies.
7. _run.py_: contains all the run parameters and variables needed for all the other .py files.
