import os
from astropy.io import fits
import multiprocessing
from joblib import Parallel, delayed
import numpy as np
import glob

from astropy.convolution import convolve


num_cores = multiprocessing.cpu_count()


def run_mfmtk(f):
    try:
        if(f in already):
            print('{} already measured'.format(f))
        else:
            path_to = base_dir + f
            print(path_to)
            data = fits.open(base_dir + f)[0]
            s = f.split('.')[0]
            os.system('python /data/morfometryka700.py {} /data/CANDELS/result00_psf.fits noshow'.format(path_to))
    except:
        print(f)
        return 0

fields = ['gs']


for field in fields:
    base_dir = '/data/CANDELS/selected_stamps/{}/'.format(field)

    files = os.listdir(base_dir)

    already_measured = glob.glob(base_dir + '*.mfmtk')
    already = [f.split('.mfmtk')[0].split('/')[-1] + '.fits'  for f in already_measured]

    print(already)
    Parallel(n_jobs=8)(delayed(run_mfmtk)(f) for f in files)


