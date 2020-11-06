import os
from astropy.io import fits
import multiprocessing
from joblib import Parallel, delayed
import numpy as np
import glob
import galckito

from astropy.convolution import convolve


num_cores = multiprocessing.cpu_count()

base_dir = '/data/TNG/mfmtk_outputs/'


files = os.listdir(base_dir)

already_measured = glob.glob('/data/TNG/mfmtk_outputs/*.mfmtk')
already = [f.split('.mfmtk')[0].split('/')[-1] + '.fits'  for f in already_measured]

print(already)

def run_mfmtk(f):
    try:
        if(f in already):
            print('{} already measured'.format(f))
        else:
            path_to = ''
            print(base_dir + f)
            data = fits.getdata(base_dir + f)
            s = f.split('.')[0]
            path_to = '/data/TNG/mfmtk_outputs/{}.fits'.format(s)
            #print(path_to)
            #fits.writeto(path_to,  data, clobber=True)
            os.system('python ~/mfmtk/morfometryka700.py {} psf.fits noshow'.format(path_to))
            os.system('rm {}'.format(path_to))
            #name = f.split('.fits')[0]

            #os.system('cat /data/mfmtk/sfr_psfed/{}_*.mfmtk >> /data/mfmtk/sfr_psfed/{}.mfmtk'.format(name, name))
            #os.system('rm /data/mfmtk/sfr_psfed/{}_*.mfmtk'.format(name))
    except (IndexError, TypeError):
        print(f)
        return 0

Parallel(n_jobs=4)(delayed(run_mfmtk)(f) for f in files)


