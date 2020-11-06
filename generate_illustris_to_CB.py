from astropy.io import fits
import pandas as pd
import  h5py
from astropy.cosmology import LambdaCDM
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import zoom

mergers = pd.read_csv('/data/captain/TNG/obs_z052/mergers.mfmtk', sep=',')
nonmergers = pd.read_csv('/data/captain/TNG/obs_z052/non_mergers.mfmtk', sep=',')


for i, row in mergers.iterrows():
    break
    src = f'/data/captain/TNG/obs_z052/cutouts/{row.rootname.strip()}.fits'
    print(row.rootname.strip(), i)
    data = fits.getdata(src)
    scaled = zoom(data, 60/128)
    try:
        fits.writeto(filename=f'/data/captain/illustris_cutouts/mergers/{row.rootname.strip()}.fits', data=scaled)
    except:
        continue

for i, row in nonmergers.iterrows():
    src = f'/data/captain/TNG/obs_z052/cutouts_nonmergers/{row.rootname.strip()}.fits'
    print(src, row.rootname.strip(), i)
    try:
       data = fits.getdata(src)
       scaled = zoom(data, 60/128)
    
       fits.writeto(filename=f'/data/captain/illustris_cutouts/non_mergers/{row.rootname.strip()}.fits', data=scaled)
    except OSError as ex:
        print(ex)
        continue
