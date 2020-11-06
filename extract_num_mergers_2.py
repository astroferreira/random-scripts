import numpy as np
import logging
import h5py
import glob
import six
import os
import sys
import argparse
import pandas as pd
from astropy.cosmology import LambdaCDM

h = 0.6774
cosmo = LambdaCDM(H0=100*h, Om0=0.3089, Ob0=0.0486, Ode0=0.6911, Tcmb0=2.73)


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


import traceback

def extract_merger_info(tree):

    df_tree = pd.DataFrame(climb_branch(tree).T, columns=['depth', 'ratios', 'origin_snap', 'central_snap', 'fp_subfind', 'np_subfind', 'fp_mass', 'np_pass'])
    return df_tree
    
def climb_branch(tree, mass_limit = 9, index=0, depth=0):
    
    rootID = tree['SubhaloID'][index]
    fpID = tree['FirstProgenitorID'][index]
    oriSnapNum = tree['SnapNum'][index]
    
    while fpID != 1:
        try:
            fpIndex = index + (fpID - rootID)
            fpMass  = maxPastMass(tree, fpIndex, 'stars')
            fpLogMass = to_log_mass(fpMass)
            fpSubfind = tree['SubfindID'][fpIndex]
            fpSnapNum = tree['SnapNum'][fpIndex]
                
            npID = tree['NextProgenitorID'][fpIndex]
            while npID != -1:

                npIndex = index + (npID - rootID)
                npMass  = maxPastMass(tree, npIndex, 'stars')
                npLogMass = to_log_mass(npMass)
                npSubfind =  tree['SubfindID'][npIndex]
                npSnapNum = tree['SnapNum'][npIndex]
                ratio = min(fpMass / npMass, npMass / fpMass)

                if ratio >= 0.25:
                    #print(depth, ratio, oriSnapNum, fpSnapNum, npSnapNum, fpSubfind, npSubfind, fpLogMass, npLogMass)
                    depths.append(depth)
                    ratios.append(ratio)
                    oriSnaps.append(str(int(oriSnapNum)))
                    mergersSnaps.append(fpSnapNum)
                    fpSubfinds.append(fpSubfind)
                    npSubfinds.append(npSubfind)
                    fpMasses.append(fpLogMass)
                    npMasses.append(npLogMass)
                
                climb_branch(tree, index=npIndex, depth=(depth+1))
                    

                npID = tree['NextProgenitorID'][npIndex]

            fpID = tree['FirstProgenitorID'][fpIndex]
        except ValueError as e:
            #print(traceback.format_exc())
            break
            
    return np.array([depths, ratios, oriSnaps, mergersSnaps, fpSubfinds, npSubfinds, fpMasses, npMasses])


def get_current_mass(tree, index, partType='stars'):

    ptNum = partTypeNum(partType)

    return tree['SubhaloMassType'][index, ptNum]

def numMergers(tree, minMassRatio=0.25, massPartType='stars', index=0):
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
    return np.round(np.log10(illustris_mass * 1e10 / cosmo.h), 3) 


import pandas as pd
import numpy as np

df = pd.read_pickle('/data/MERGERS/datasets/df_sample_with_sfr_compact.pk')
#only_mergers = df.iloc[np.where(df.merger_label < 2)]
#z0_subfind = only_mergers.z0_subfind

dfs = []
for j, (idx, row) in enumerate(df.iterrows()):
    if j % 10 == 0:
        print(j / df.shape[0])

    depths = []
    ratios = []
    oriSnaps = []
    mergersSnaps = []
    fpSubfinds = []
    npSubfinds = []
    fpMasses = []
    npMasses = []
    tree = h5py.File(f'/data/MERGERS/merger_trees/{str(int(row.snap))}_{str(int(row.subfind))}.hdf5', 'r')
    df_tree = extract_merger_info(tree)
    df_tree['snap'] = str(int(row.snap))
    df_tree['subfind'] = str(int(row.subfind))
    df_tree['rootname'] = f'{str(int(row.snap))}_{str(int(row.subfind))}'
    #print(df_tree)
    dfs.append(df_tree)

df_complete = pd.concat(dfs)
print(df_complete)
df_complete.to_pickle('merger_info_2.pk')