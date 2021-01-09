Dependencies:
1. python 3+ (anaconda recommended; conda create -n nupyprop conda)
2. numpy (conda install numpy)
3. pandas (conda install pandas)
4. scipy (conda install scipy)
5. matplotlib (conda install matplotlib)
6. numba (conda install numba)
7. sympy (conda install sympy)
8. interpolation (conda install -c conda-forge interpolation)
9. (optional for GUI HDF5 viewer) vitables (conda install -c conda-forge vitables)

Run parameters are defined in run.py. Different switches are described as follows:

1. -e or --energy: incoming neutrino energy in GeV. Works for single energy or multiple energies. For multiple energies, separate energies with commas eg. 1e7,1e8,1e9. Default energies are 1e7,1e8,1e9,1e10,1e11 GeV.

2. -a or --angle: slant Earth angles in degrees. Works for single angle or multiple angles. For multiple angles, separate angles with commas eg. 1,3,5,7,10. Default angles are 1,3,5,7,10,12,15,17,20,25,30,35 degrees.

3. -i or --idepth: depth of ice/water in km. Default value is 4 km.

4. -l or --lepton: flavor of lepton used to propagate. Can be either muon or tau. Default is tau.

5. -m or --material: material used in electromagnetic energy loss; not used in main program, only used for running energy_loss.py individually. Default is rock.

6. -x or --xc_model: neutrino/anti-neutrino cross-section model used. For a list of model names, see lookup_tables.h5/Neutrino_Cross_Sections. Default is ncteq15. 

7. -p or --pn_model: photonuclear energy loss model used. For now, can be either bb (Bezrukov-Bugaev) or allm (Abramowicz, Levin, Levy, Maor). Default is allm.

8. -f or --fac_nu: rescaling factor for SM cross-sections. Default is 1.

9. -s or --stat: statistics. Default is 1e7 particles.

Example for running tau propagation for 10^8 GeV neutrinos at 10 degrees with a statistics of 10^8 particles with all other parameters as defaults:
python run.py -e 1e8 -a 10 -s 1e8

Viewing output results:
output.h5 will contain the results of the code after it has finished running. 

Currently, the output file contains a full set of tau results ran for neutrino energies of 1e7, 1e8, 1e9 and 1e10 GeV, at 1,3,5,7,10,12,15,17,20,25,30,35 degrees with statistics of 1e8 particles (10^9 GeV ran with 10^7 statistics because of impending memory leak errors.)

Note: Since the HDF5 files have been written using pandas, the full tables are stored in 'axis1' in all of the branches.

WARNING: Running the code will replace the results of the output.h5 file. Backing up previous output files is recommended (or just download and replace from the repo). Future fixes include overwrite warnings for the user.

For debugging/development:
The correct order to look at the code is in the following order:

1. data.py: contains the reading/writing modules from/to the hdf5 files.
2. geometry.py: contains the Earth geometry modules (including PREM).
3. cross_section.py: contains neutrino/anti-neutrino cross_section models.
4. energy_loss.py: contains lepton energy loss models.
5. my_interpolation.py: contains (numba enabled) interpolation functions.
6. transport.py: heart of the code; contains functions to propagate neutrinos and leptons through the Earth.
7. main.py: forms the main skeleton of the code; propagates the neutrinos and leptons, and calculates the p_exit and collects outgoing lepton energies.
8. run.py: contains all the run parameters and variables needed for all the other .py files.