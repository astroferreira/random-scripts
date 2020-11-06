from astropy.cosmology import LambdaCDM
import pandas as pd
import numpy as np
h = 0.6774
cosmo = LambdaCDM(H0=100*h, Om0=0.3089, Ob0=0.0486, Ode0=0.6911, Tcmb0=2.73)


import os

"""

Illustris Public Script helper functions

"""
def maxPastMass(tree, index, partType='stars'):
    """ Get maximum past mass (of the given partType) along the main branch of a subhalo
        specified by index within this tree. """
    ptNum = partTypeNum(partType)

    branchSize = tree['MainLeafProgenitorID'][index] - tree['SubhaloID'][index] + 1
    masses = tree['SubhaloMassType'][index: index + branchSize, ptNum]
    return np.max(masses)

def partTypeNum(partType):
    """ Mapping between common names and numeric particle types. """
    if str(partType).isdigit():
        return int(partType)
        
    if str(partType).lower() in ['gas','cells']:
        return 0
    if str(partType).lower() in ['dm','darkmatter']:
        return 1
    if str(partType).lower() in ['tracer','tracers','tracermc','trmc']:
        return 3
    if str(partType).lower() in ['star','stars','stellar']:
        return 4 # only those with GFM_StellarFormationTime>0
    if str(partType).lower() in ['wind']:
        return 4 # only those with GFM_StellarFormationTime<0
    if str(partType).lower() in ['bh','bhs','blackhole','blackholes']:
        return 5
    
    raise Exception("Unknown particle type name.")

cat = np.load('/data/snapTNG.npy')
gys = cat[2]
snapshots= cat[0]





def get_current_mass(tree, index, partType='stars'):

    ptNum = partTypeNum(partType)

    return tree['SubhaloMassType'][index, ptNum]

def numMergers(tree, minMassRatio=1e-01, massPartType='stars', index=0):
    """ Calculate the number of mergers in this sub-tree (optionally above some mass ratio threshold). """
    # verify the input sub-tree has the required fields
    reqFields = ['SubhaloID', 'NextProgenitorID', 'MainLeafProgenitorID',
                 'FirstProgenitorID', 'SubhaloMassType']

    central_snapshots = np.array([])
    central_subhalos = np.array([])

    merger_next_progenitors = np.array([])
    merger_first_progenitors = np.array([])
    merger_ratios = np.array([])
    fp_log_masses = np.array([])
    np_log_masses = np.array([])

    if not set(reqFields).issubset(tree.keys()):
        raise Exception('Error: Input tree needs to have loaded fields: '+', '.join(reqFields))

    numMergers   = 0
    invMassRatio = 1.0 / minMassRatio

    # walk back main progenitor branch
    rootID = tree['SubhaloID'][index]
    fpID   = tree['FirstProgenitorID'][index]

    current_central_snapshot = tree['SnapNum'][index]
    current_central_subhalo = tree['SubfindID'][index]

    while fpID != -1:

        fpIndex = index + (fpID - rootID)
        fpMass  = maxPastMass(tree, fpIndex, massPartType)

        # explore breadth
        npID = tree['NextProgenitorID'][fpIndex]

        while npID != -1:
            npIndex = index + (npID - rootID)
            npMass  = maxPastMass(tree, npIndex, massPartType)
            
            np_log_mass = to_log_mass(npMass)
            fp_log_mass = to_log_mass(fpMass)

            mass_limit = 8.0

            if fp_log_mass > mass_limit or np_log_mass > mass_limit:
            
                ratio = npMass / fpMass
                
                if ratio >= minMassRatio and ratio <= invMassRatio:
                    numMergers += 1

                    central_snapshots = np.append(central_snapshots, current_central_snapshot)
                    central_subhalos = np.append(central_subhalos, current_central_subhalo)

                    merger_next_progenitors = np.append(merger_next_progenitors, npIndex)
                    merger_first_progenitors= np.append(merger_first_progenitors, fpIndex)
                    merger_ratios = np.append(merger_ratios, ratio)
                    fp_log_masses = np.append(fp_log_masses, fp_log_mass)
                    np_log_masses = np.append(np_log_masses, np_log_mass)

            npID = tree['NextProgenitorID'][npIndex]
	
	
        fpID = tree['FirstProgenitorID'][fpIndex]
        current_central_subhalo = tree['SubfindID'][fpIndex]
        current_central_snapshot = tree['SnapNum'][fpIndex]
        
    return (numMergers, merger_first_progenitors.astype(int), merger_ratios, fp_log_masses, np_log_masses, merger_next_progenitors.astype(int), central_snapshots, central_subhalos)

def to_log_mass(illustris_mass):
    return np.round(np.log10(illustris_mass * 1e10 / cosmo.h), 2) 


snapZ = np.load('/data/snapTNG.npy')
headers = {"api-key":"ff352a2affacf64753689dd603b5b44e"}
import requests
def get(path, filename=None, params=None):
    # make HTTP GET request to path

    if(~os.path.exists(f'/data/MERGERS/merger_trees/{filename}.hdf5')):
        r = requests.get(path, params=params, headers=headers)

        # raise exception if response code is not HTTP SUCCESS (200)
        r.raise_for_status()

        if r.headers['content-type'] == 'application/json':
            return r.json() # parse json responses automatically
        
        if 'content-disposition' in r.headers:
            #print(r.headers['content-disposition'].split("filename=")[1])
            with open('/data/MERGERS/merger_trees/'+filename+'.hdf5', 'wb') as f:
                f.write(r.content)
                f.flush()
                f.close()
            return '/data/MERGERS/merger_trees/'+filename+'.hdf5' # return the filename string
    else:
        return 0

    return r

import sys

df = pd.read_pickle('/data/MERGERS/datasets/df_sample_with_sfr_compact.pk')

df = df.iloc[np.where(df.merger_label == int(sys.argv[1]))]


for i, row in df.iterrows():
    #print(i)
    path = f'http://www.tng-project.org/api/TNG300-1/snapshots/{str(int(row.snap))}/subhalos/{str(int(row.subfind))}/sublink/full.hdf5'    
    filename = f'{str(int(row.snap))}_{str(int(row.subfind))}'
    if(os.path.exists(f'/data/MERGERS/merger_trees/{filename}.hdf5')):
        #print('skip')
        continue
    try:
        get(path, filename=filename)
    except Exception as e:
        print(e)
        continue