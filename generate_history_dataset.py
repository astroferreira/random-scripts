import pandas as pd
import numpy as np

df = pd.read_pickle('/data/MERGERS/datasets/df_sample_with_sfr.pk')
only_mergers = df.iloc[np.where(df.merger_label < 2)]
z0_subfind = only_mergers.z0_subfind

unique_at_z0 = np.unique(z0_subfind)
n_merging_events = np.zeros_like(unique_at_z0)

for i, subfind in enumerate(unique_at_z0):
    temp = only_mergers.iloc[np.where(z0_subfind == subfind)]
    n_merging_events[i] = temp.shape[0]

unique, counts = np.unique(n_merging_events, return_counts=True)
pd.DataFrame(np.array([unique, counts]).T, columns=['MM in the Past', 'Counts'])

columns=['snap', 'subfind', 'np_subfind', 'z0_subfind', 'ratio', 'fp_mass', 'np_mass', 'central_snap', 'central_subfind', 'num', 'matches']
tng100_mergers = pd.read_csv('/data/captain/sources/tireless-tracker/TNG100_mergers.dat', sep=';', names=columns)


z0_subfinds = np.unique(tng100_mergers.z0_subfind)
num_mergers = np.zeros_like(z0_subfinds)
for i, sub in enumerate(z0_subfinds):
    num_mergers[i] = tng100_mergers.iloc[np.where((tng100_mergers.z0_subfind == sub) & (tng100_mergers.ratio > 0.25))].shape[0]


available = np.loadtxt('/data/captain/TNG/sdss/snapnum_099/subfind_ids.txt').astype(int)



reference = pd.DataFrame(np.array([z0_subfinds, num_mergers]).T, columns=['subfind', 'Counts'])
reference = reference.set_index('subfind')
availables = reference.reindex(available)
availables.iloc[np.where(np.isnan(availables.Counts))] = 0
print(availables.groupby('Counts').size())


import glob
broadbands = glob.glob('/data/captain/TNG/sdss/snapnum_099/data/*.fits')


from IPython.display import Image
from astropy.io import fits 
from matplotlib import pyplot as plt
from scipy.ndimage import zoom

dataset = np.zeros((availables.shape[0], 128, 128, 4))

for k, (index, row) in enumerate(availables.iterrows()):
    
    if(k % 100 == 0):
        print(k)
        
    data = fits.getdata(f'/data/captain/TNG/sdss/snapnum_099/data/broadband_{index}.fits')

    g = data[0]
    g = zoom(g, 128/data[0].shape[0])
    r = data[1]
    r = zoom(r, 128/data[1].shape[0])
    i = data[2]
    i = zoom(i, 128/data[2].shape[0])
    z = data[3]
    z = zoom(z, 128/data[3].shape[0])
    
    dataset[k] = np.dstack([g, r, i, z])

np.save(arr=dataset, file='dataset.npy')