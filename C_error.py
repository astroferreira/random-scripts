import scipy
import pandas as pd
import numpy as np

from astropy.io import fits
from matplotlib import pyplot as plt
from photutils import aperture_photometry, CircularAperture

def _Concentration(R20, R80):
    return 5*np.log10(R80/R20)

def estimate_C_error(galaxy, row, exp_time):

    R10=CircularAperture([row.x0peak, row.y0peak], row.R10)
    R20=CircularAperture([row.x0peak, row.y0peak], row.R20)
    R30=CircularAperture([row.x0peak, row.y0peak], row.R30)
    R40=CircularAperture([row.x0peak, row.y0peak], row.R40)
    R50=CircularAperture([row.x0peak, row.y0peak], row.R50)
    R60=CircularAperture([row.x0peak, row.y0peak], row.R60)
    R70=CircularAperture([row.x0peak, row.y0peak], row.R70)
    R80=CircularAperture([row.x0peak, row.y0peak], row.R80)
    R90=CircularAperture([row.x0peak, row.y0peak], row.R90)
   
    fracs = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    radii = np.array([R10, R20, R30, R40, R50, R60, R70, R80, R90])
    radiipix = np.array([row.R10, row.R20, row.R30, row.R40, row.R50, row.R60, row.R70, row.R80, row.R90])

    table = aperture_photometry(galaxy, radii)

    t = np.array(table)[0]
    LTs = np.array([t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10], t[11]])*exp_time

    shot20 = np.sqrt(LTs[1])
    shot80 = np.sqrt(LTs[7])


    LR20 = scipy.interpolate.splrep(radiipix, LTs - (LTs[1]+shot20), s=0)
    LR80 = scipy.interpolate.splrep(radiipix, LTs - (LTs[7]-shot80), s=0)
    RR20 = scipy.interpolate.sproot(LR20)[0]     
    RR80 = scipy.interpolate.sproot(LR80)[0]   

    C1 = _Concentration(row.R20, row.R80)
    Cnoise = _Concentration(RR20, RR80)

    return C1-Cnoise  


if __name__ == '__main__':    
"""
    Example usage with MFMTK data in a Pandas DataFrame
"""

dataset = pd.read_pickle('dataset.pk')
BASE_PATH = '/data'

Cerrors = np.zeros_like(dataset.C1)
IDs = np.zeros_like(dataset.C1)

for idx, row in GN.iterrows():
    galaxy = fits.getdata(f'{BASE_PATH}/{row.rootname.strip()}.fits')
    Cerrors[idx] = estimate_C_error(galaxy, row, 1200)
    
dataset['Cerror'] = Cerrors