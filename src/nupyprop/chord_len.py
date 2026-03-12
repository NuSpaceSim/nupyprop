import numpy as np
import matplotlib.pyplot as plt

R_earth = 6371
beta_deg = 3
x = np.arange(0.1, 666, 1)
idepth = 4

#Earth layers considered in PREM - Earth density model
Rlay = np.array([1221.5, 3480.0, 5701.0, 5771.0, 5971.0, 6151.0, 6346.6, 6356.0, 6368.0, 6371.0])

def premdensity(Rin, idepth):
    '''
    PREM Earth density model. Computes density at multiple radii of Earth simultaneously. 
    hep-ph/9512364v1 eq. 25

    Parameters
    ----------
    Rin : float
        Radius (in km) at which density is to be calculated.
    idepth : integer
        Depth (in km) of water layer (sets the last layer).

    Returns
    -------
    density : float
        Density (in g/cm^3) at the specified radius
    '''
    # Update last layer for water depth
    Rlay[8] = 6368.0 + (3.0 - int(idepth))
    y = Rin / R_earth  # Normalize by Earth radius

    # Define density expressions as NumPy functions
    densities = np.array([
        13.0885 - 8.8381 * y**2,
        12.5815 - 1.2638 * y - 3.6426 * y**2 - 5.5281 * y**3,
        7.9565 - 6.4761 * y + 5.5283 * y**2 - 3.0807 * y**3,
        5.3197 - 1.4836 * y,
        11.2494 - 8.0298 * y,
        7.1089 - 3.8045 * y,
        2.6910 + 0.6924 * y,
        np.full_like(y, 2.900),
        np.full_like(y, 2.600),
        np.full_like(y, 1.020 if idepth > 0 else 2.600)
    ])

    # Use np.select() to apply the correct density function based on conditions
    conditions = [(Rin <= Rlay[i]) | ((i == len(Rlay) - 1) & (Rin <= Rlay[-1] * 1.001)) for i in range(len(Rlay))]
    edens = np.select(conditions, densities, default=0.0)  # Default to 0 if no condition matches

    return edens

chord_length = 2*R_earth*np.sin(beta_deg*(np.pi/180)) # 2 R_E sin(beta)
print("chord_len = ", chord_length)
r2 = (chord_length-x)**2 + R_earth**2 - 2*R_earth*(chord_length-x)*np.sin(beta_deg*(np.pi/180)) # just the law of cosines

if (beta_deg < 5):
    r = R_earth*(1 + 0.5*(x**2-chord_length*x)/R_earth**2)
else:
    r = np.sqrt(r2)

print("r = ", r)
dens = premdensity(r, idepth)

plt.plot(x, r)
plt.show()
plt.plot(r, dens, ls='None', marker='.')
plt.show()
plt.plot(x, dens)
plt.show()

#print(premdensity(r, idepth))