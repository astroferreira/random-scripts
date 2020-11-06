import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy import units as u
from matplotlib import pyplot as plt
from astropy.wcs.wcs import NoConvergence
from astropy.nddata.utils import NoOverlapError
from PIL import Image

import sys

df_path = sys.argv[1]
field_path = sys.argv[2]
filter = sys.argv[3]
prefix = sys.argv[4]
field_name = sys.argv[5]

print(sys.argv)

df = pd.read_pickle(df_path)
coords = [SkyCoord(ra=row.RA * u.deg, dec=row.DEC * u.deg) for i, row in df.iterrows()]

field = fits.open(field_path, memmap=True)

try:
    field[0].header.pop('CPDIS1')
    field[0].header.pop('CPDIS2')
    field[0].header.pop('D2IMDIS1')

    field[0].header.pop('DP1')
    field[0].header.pop('DP2')
except Exception as e:
    print(e)

wcs = WCS(field[0].header)
wcs.sip = None

for (i, row), coord in zip(df.iterrows(), coords):
    print(i)
    try:
        cutout = Cutout2D(field[0].data, position=coord, size=101, wcs=wcs)
    except NoOverlapError:
        pass
    id = int(row.ID)
    fits.writeto(data=cutout.data, filename=f'/data/CANDELS/cutouts/{field_name}/{filter}/fi/{prefix}_{id}.fits', overwrite=True)
