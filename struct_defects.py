#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: wwwennie

  see pymatgen documentation for more information
 
  distort the structure slightly to bias the calculation 
  towards polaron formation
  i.e., expanding/shrinking the nn bonds for electron/hole polarons

----------------------------------------------
"""

import math
from itertools import combinations
import numpy as np

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

#===============================================================================#


#------------------------- bias structure to polaron formation  --------------------------#

def bias_polaron(struc,site,maxr,bias=1.01):
   """ for a particular site in Pymatgen structure, bias the nearest neighbors w/i distance
       maxr by bias

       struct (Pymatgen structure)
       site   (integer): index of atomic site
       maxr   (float): max distance to search for nn for atom site (Angstroms)
       bias   (float): how much to change nn distances by, multiplicative factor
                           default to slight expansion for electron polaron formation
                           0.0 = no change

   """  
   # get nearest neighbors
   psite=struc.sites[site-1] # zero-indexing in python
   nn = struc.get_neighbors(psite,maxr,include_index=True)
   bias_struct = struc.copy()

   print(psite.specie,site,psite.frac_coords)

   # find new distances for neighbors
   for inn in nn:
      nnsite,nndist,nnind = inn
      
      # calculate difference vector
      diff = np.subtract(nnsite.coords,psite.coords)
      pdiff = diff * bias
      newcoord = np.add(psite.coords,pdiff) 

      # replace with new coord
      bias_struct.replace(nnind,nnsite.specie,coords=newcoord,coords_are_cartesian=True)
      
   # return new structure
   return bias_struct
 
def vacancy(structure,sites):
   """ wrapper for making vacancies 
 
       struc (Pymatgen structure)
       sites (list of integers): sites to remove """
   # copy structure and leave original unmodified
   struc = structure.copy()



   sites = np.subtract(sites,1) # zero-indexing
   struc.remove_sites(sites)
   
   return struc


def substitute(structure,sites,species):
   """ wrapper for making substitutional defects
 
       struc (Pymatgen structure)
       sites (list of integers): sites to replace 
       species (list of strings): species of sites """
   # copy structure and leave original unmodified
   struc = structure.copy()


   sites = np.subtract(sites,1) # zero-indexing
   assert len(sites) == len(species)

   for i in range(len(sites)):
      s = sites[i]
      coord = struc.sites[s].frac_coords
      struc.replace(s,species[i],coord)
  
   return struc

def translate(structure,sites,vec_trans):
   """ translate site by vec_trans

       hackneyed attempt to make the structure more "corner sharing"
       see Seo2018 Role of Point Defects
 
       vec_trans (list of floats): containing delta_x, delta_y, delta_z; 
                                   assume fractional coord
       sites (list of integers)

   """
   # copy structure and leave original unmodified
   struc = structure.copy()

   for site in sites:
 
     specie = struc.sites[site].specie
     atom = struc.sites[site].frac_coords

     newcoord = np.add(atom,vec_trans)
     struc.replace(site,specie,newcoord)

   return struc 

def translate_to(structure,site,targetsite,shift=0.15):
   """ translate site by to targetsite by multiplicative amount shift

       hackneyed attempt to make the structure more "corner sharing"
       see Seo2018 Role of Point Defects
       sites, targetsite (list of integers, integer): corresponding to site

   """
   # copy structure and leave original unmodified
   struc = structure.copy()

 
   site = site - 1  # zero-indexing
   targetsite = targetsite - 1

   specie = struc.sites[site].specie
   atom = np.mod(struc.sites[site].frac_coords,1) # periodic boundary conditions
   target = np.mod(struc.sites[targetsite].frac_coords,1)
   vec_trans = np.multiply(np.subtract(target,atom),shift)

   newcoord = np.add(atom,vec_trans)
   struc.replace(site,specie,newcoord)

   return struc 

def translate_nn(struc,site,maxr,target,shift=0.3):
   """ Translate atom indexed by site in sites and nn within maxr dist to coordinates specified
        by target index
  
       sites(list of integers): corresponding to sites, automatic zero-index calibrated
       target(integer): corresponding to site with coordinates, e.g., towards VO3 + V_O       
       maxr (float): radial distance to search to nn
       shift (float): amount to shift by, multiplicative factor of dist to target

       yet another attempt to stabilize polaron"""

   newstruc = struc.copy()

   newstruc = translate_to(newstruc,site,target,shift)

   nn = struc.get_neighbors(struc.sites[site],maxr,include_index=True)
   for inn in nn:
      nnsite,nndist,nnind = inn
      newstruc = translate_to(newstruc,nnind+1,target,shift)

   return newstruc

def translate_nn_rigid(struc,site,maxr,target,shift=0.3):
   """ Translate atom indexed by site in sites and nn within maxr dist to coordinates specified
        by target index
       Rigid shift of atom + nn
  
       site(integer): corresponding to sites, automatic zero-index calibrated
       target(integer): corresponding to site with coordinates, e.g., towards VO3 + V_O       
       maxr (float): radial distance to search to nn
       shift (float): amount to shift by, multiplicative factor of dist to target

       yet another attempt to stabilize polaron"""

   newstruc = struc.copy()
   site = site - 1
   target = target - 1
 
   coord = struc.sites[site].frac_coords
   targetcoord = struc.sites[target].frac_coords
   vec_trans = np.multiply(np.subtract(targetcoord,coord),shift)
   #print(targetcoord, coord)
   #print(vec_trans)

   nn = struc.get_neighbors(struc.sites[site],maxr,include_index=True)
   newstruc = translate(newstruc,[site],vec_trans)

   for inn in nn:
      nnsite,nndist,nnind = inn
      newstruc = translate(newstruc,[nnind],vec_trans)

   return newstruc

def rotate(structure,sites,angle):
   """ Rotate site1 and site2 about average point by angle 
          in xy-plane 
  
       sites(list of integers): corresponding to sites, automatic zero-index calibrated
       angle (degrees)        

       yet another attempt to stabilize polaron"""

   # copy structure and leave original unmodified
   struc = structure.copy()

   # at the current moment, only need rotation within xy-plane
   R = lambda ang: np.asmatrix([[np.cos(ang),-np.sin(ang), 0],[np.sin(ang),np.cos(ang),0],[0,0,1]])
   sites = np.subtract(sites,1) # zero-indexing

   # simple average coordinate; untested for lots of atoms of different masses
   avg = np.array([0,0,0])
   for s in sites:      
      site = struc.sites[s].frac_coords
      avg = np.add(avg,site)
   avgsite = np.divide(avg,len(sites))

   # rotation about arbitrary center
   # T(x,y) R T(-x,-y)
   for s in sites:
      # rotate coordinates
      site = struc.sites[s].frac_coords
      specie = struc.sites[s].specie
      newcoord = R(np.radians(angle)) * np.subtract(site,avgsite).reshape((3,1))
      newcoord = np.add(newcoord.reshape((1,3)),avgsite)
      #print(newcoord.A1)
   
      # substitute coordinates
      struc.replace(s,specie,newcoord.A1)
  
   return struc

