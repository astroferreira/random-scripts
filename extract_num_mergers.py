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
    # verify t