Dependencies:
1. python 3+ (anaconda recommended; conda create --n nupyprop conda)
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