#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: wwwennie

  see pymatgen documentation for more information
 
  self-contained module of helper functions for 
  simulating STM images 
----------------------------------------------
"""

import math
from itertools import combinations
import numpy as np

from struct_import import *
#------------------------------- Helper functions --------------------------------#

def cryst2cart(vec,trmat,iflag):

  """
    convert from crystalline to cartesian coordinates and back

    vec = vector
    trmat = transformation matrix
    if iflag = 1, crystallographic to cartesian
        trmat = real-space lattice basis for atoms
              = reciprocal-space lattice basis for k-point
    if iflag = -1: the opposite 
  """
  
  vau = np.zeros(3)

  # iflag cases got flipped for single vector (compared to original in QE)
  # likely since QE version does implicit transpose for the case of many vectors
  #   checked explicitly, but didn't think too much about why 
  if iflag == -1:
     for kpol in range(3):
        vau[kpol] = trmat[kpol][0] * vec[0] + trmat[kpol][1]*vec[1] + trmat[kpol][2]*vec[2]
  elif iflag == 1:
     for kpol in range(3):
        vau[kpol] = trmat[0][kpol] * vec[0] + trmat[1][kpol]*vec[1] + trmat[2][kpol]*vec[2]
  else:
     print("I need a valid iflag option!")

  return vau

def getCartCoord(struc,natoms):
    """Extracting x,y,z Cartesian coordinates in pymatgen Structure"""
    # struc = structure from pymatgen
    # natoms = number of atoms
    x=[]
    y=[]
    z=[]
    for nat in range(natoms):
        x.append(struc.sites[nat].x)
        y.append(struc.sites[nat].y)
        z.append(struc.sites[nat].z)
    return x,y,z

def get_surf_atoms(filename,spec,nlayer=2,is_symm_slab=True):
   """
     Wrapper for id_surf_atoms to handle QE inputs and give output for bash scripts to parse
     Search for surf atoms of species spec within nlayers    

     filename: (string) QE input file
     spec: (string) species element
     nlayers: (integer) number of layers to search through
     is_symmm_slab: (boolean) if symmetric cell, return surface atoms on both exposed surfaces
                      otherwise search for atoms within max +z coordinates

     surfatoms: list of atom indices corresponding to surface (python zero-based indexing)
 
     as with struct_import, depends on pymatgen and installation of XCrysden
     more specifically on pwi2xsf and pwo2xsf (called upon from qe2pmg.sh)
   """
   
   struc = struct_import(filename) 
   # list of Atom objects and their indices
   surf_Atoms,surfatoms_ind = id_surf_atoms(struc,spec,nlayer,is_symm_slab)
   
   return surf_Atoms, surfatoms_ind


#===============================================================================#

#-------------------  @ surface ------------------------#
def id_surf_atoms(struc,spec1,n=2,is_symm_cell=True):
    """ 
    look for atoms within n layers of the surface of species spec and return their indices
 
    Input
     struc = structure object from pymatgen
     natoms = number of atoms, integer
     maxr = radius for searching distances for each site i, float;units of Ang.
     spec1= identity of speciesi,  chemical symbol; string
     n = integer, number of layers from surface to consider; +/- z
     is_symm_cell = if slab is symmetric, returns surface atoms on both exposed surfaces
    Output
     surfatoms: list of Atoms
    """
    
    natoms = struc.num_sites

    # for filtering distances from only species 1
    speclist=np.array(struc.indices_from_symbol(spec1))

    #assumes slab constructed long-way along z
    #find z-cutoff for atoms of species 1 considered at surface 
    zcoords = np.zeros(len(speclist))
    for ispec in np.arange(len(speclist)):
      spec = speclist[ispec]
      zcoords[ispec] = struc.sites[spec].coords[2]

    # +/- z, average of n and n+1 inner layer z coordinates 
    uzcoords = np.unique(zcoords.round(decimals=1)) # implicit tolerance
    uzcoords.sort()
    mz = np.average([uzcoords[n],uzcoords[n-1]])
    pz = np.average([uzcoords[-n],uzcoords[-n-1]])
    print("searching z <", mz)
    print("searching z >",pz)

    surfatoms = []
    surfatoms_ind = []
    for icat in range(natoms): 
      # check within +/- z of surface  
      z = struc.sites[icat].coords[2]
      if (z < mz):
          # check for species
          if (icat in speclist):
             surfatoms.append(struc.sites[icat])
             surfatoms_ind.append(icat)
      if is_symm_cell and (z > pz):
          if (icat in speclist):
             surfatoms.append(struc.sites[icat])
             surfatoms_ind.append(icat)
    return surfatoms, surfatoms_ind
 
