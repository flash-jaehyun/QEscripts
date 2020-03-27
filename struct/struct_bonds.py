#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: wwwennie

  see pymatgen documentation for more information
 
  find bond lengths and bond angles
  species-dependent selection, e.g., look for only A-B-A kinds of bonds
  periodic boundary conditions included

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


#------------------------- Single-instance functions --------------------------#

def getlength(struc,maxr):
    """Finding number nearest neighbors of distance r    
      search for atoms within r radius of each atomic site

     Input
      struc = structure read from pymatgen
      maxr = radius for searching distances for each site i, float; units of Ang.
     Output
      rad = list of distances within maxr radius
    """
   
    natoms = struc.num_sites  

    # leads to double counting of bonds?
    rad = []
    neigh_list=struc.get_all_neighbors(r=maxr,include_index=False)
    for icat in range(natoms): 
        #print len(neigh_list[itr])
        n_num = len(neigh_list[icat])
        for ineigh in range(n_num):
            dist = neigh_list[icat][ineigh][1]
            rad.append(dist);
    return rad

def getangle_sites(struc,site1,site2,site3):
    """   
      get angle between sites  

      included is when site is near periodic boundary edge; searches for image of 
      neighbor
 
      pymatgen also has get_angle(), but it does not appear to include boundary 
        conditions, even though it is implemented elsewhere
    
      Input:
        structure = pymatgen structure
        site1, site2, site3 = pymatgen Sites; assumes site2 is center atom


      differs from get_angle in that inputs are Sites
    """

    v1 = site1.coords-site2.coords
    v2 = site3.coords-site2.coords

    v1cryst = site1.frac_coords-site2.frac_coords
    v2cryst = site3.frac_coords-site2.frac_coords

    # for the purposes of applying periodic boundary conditions
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
  
    # periodic boundary conditions
    #  smallest norm is assumed to be the proper periodic image of the site
    translate = np.identity(3)
    for i in range(3):
      # uses fractional coordinates, but returns in absolute distance
      # (tmpnorm1,trans_vec1) = site1.distance_and_image(site2,jimage=translate[i])       
      # (tmpnorm2,trans_vec2) = site3.distance_and_image(site2,jimage=translate[i])       
      
      (tmpnorm1,trans_vec1) = site1.distance_and_image(site2)       
      (tmpnorm2,trans_vec2) = site3.distance_and_image(site2)       

      if (tmpnorm1  < norm1):
         norm1 = tmpnorm1
         # move site1 to be within same cell as site2
         v1cryst = (site1.frac_coords - trans_vec1) - site2.frac_coords 
         v1 = cryst2cart(v1cryst,struc.lattice.matrix,1)
     
      if (tmpnorm2 < norm2):
         norm2 = tmpnorm2
         # move site3 to be within same cell as site2
         v2cryst = (site3.frac_coords - trans_vec2) - site2.frac_coords
         v2 = cryst2cart(v2cryst,struc.lattice.matrix,1)
     
    # calculate angle, from pymatgen.util.coord 
    d = np.inner(v1, v2) / norm1 / norm2
    d = min(d, 1)
    d = max(d, -1)
    angle = math.acos(d)

    return math.degrees(angle) 

def getangle(struc,i1,i2,i3):
    """   
      get angle between sites i1, i2 (middle), i3 (python indexing) 

      included is when site is near periodic boundary edge; searches for image of 
      neighbor
 
      pymatgen also has get_angle(), but it does not appear to include boundary 
        conditions, even though it is implemented elsewhere
    
      Input:
        structure = pymatgen structure
        i1, i2, i3 = atomic coordinate indices; i2 is assumed center atom

     differs from getangle_sites in that inputs are by site indices of struc
     
    """

    site1=struc.sites[i1]
    site2=struc.sites[i2]
    site3=struc.sites[i3]
   
    v1 = site1.coords-site2.coords
    v2 = site3.coords-site2.coords

    v1cryst = site1.frac_coords-site2.frac_coords
    v2cryst = site3.frac_coords-site2.frac_coords

    # for the purposes of applying periodic boundary conditions
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
  
    # periodic boundary conditions
    #  smallest norm is assumed to be the proper periodic image of the site
    translate = np.identity(3)
    for i in range(3):
      # uses fractional coordinates, but returns in absolute distance
      # (tmpnorm1,trans_vec1) = site1.distance_and_image(site2,jimage=translate[i])       
      # (tmpnorm2,trans_vec2) = site3.distance_and_image(site2,jimage=translate[i])       
      
      (tmpnorm1,trans_vec1) = site1.distance_and_image(site2)       
      (tmpnorm2,trans_vec2) = site3.distance_and_image(site2)       

      if (tmpnorm1  < norm1):
         norm1 = tmpnorm1
         # move site1 to be within same cell as site2
         v1cryst = (site1.frac_coords - trans_vec1) - site2.frac_coords 
         v1 = cryst2cart(v1cryst,struc.lattice.matrix,1)
     
      if (tmpnorm2 < norm2):
         norm2 = tmpnorm2
         # move site3 to be within same cell as site2
         v2cryst = (site3.frac_coords - trans_vec2) - site2.frac_coords
         v2 = cryst2cart(v2cryst,struc.lattice.matrix,1)
     
    # calculate angle, from pymatgen.util.coord 
    d = np.inner(v1, v2) / norm1 / norm2
    d = min(d, 1)
    d = max(d, -1)
    angle = math.acos(d)

    return math.degrees(angle) 

#====================================================================================#

#------------------------ Species specific distribution -----------------------------#
def spec_bond_len(struc,maxr,spec1,spec2):
    """Finding number of nearest Neighbors of distance r, version 2, with species distinction
    search for atoms within r radius (I think in units of lattice vectors) of each atomic site
    and bin the distances in histogram

    get unique number of bond lengths
    give maxr slightly more than 2nd nn distance, 
    identify by element species 
 
    Input
     struc = structure read from pymatgen
     natoms = number of atoms, integer
     maxr = radius for searching distances for each site i, float;units of Ang.
     spec1, spec2 = identity of species associated with distance, chemical symbol; string
    Output
     rad: list of bond distances between spec1 and spec2

    """
    
    rad = []
    natoms = struc.num_sites
    
    # get distances of neighbors within distance maxr
    neigh_list=struc.get_all_neighbors(r=maxr,include_index=False)
    # for filtering distances from only species 1
    speclist=struc.indices_from_symbol(spec1)

    for icat in range(natoms): 
        # check for species
        if (icat in speclist):
            #print len(neigh_list[itr])
            n_num = len(neigh_list[icat])
            for ineigh in range(n_num):
                point=neigh_list[icat][ineigh]
                #print point
                #check if finding distance to spec 2 
                tmpspec2=point[0].as_dict()['species'][0]['element']
                if (tmpspec2 == spec2):
                    dist = point[1]
                    rad.append(dist);
                    #print spec1,spec2,tmpspec2,dist
    return rad
 
def spec_bond_angle(struc,maxr,spec1,spec2,minangle=0.0,maxangle=180.0):
    """Finding bond angle distribution with spec 1 as the center atom, crystalline with definite 
        coordination number
    search for atoms within r radius (I think in units of lattice vectors) of each atomic site
    and bin the distances in histogram

    Inputs
     struc = structure read from pymatgen
     maxr = radius for searching distances for each site i, float; I think in units of Ang.
     spec1, spec2 = identity of species associated with distance; string
     minangle, maxangle = bounds of angles looking for
    Outputs
     ang = list of angles between for bonds of spec2-spec1-spec2
    """
    
    #from . import getangle # needed? 
  
    ang = []
    natoms = struc.num_sites
    
    # get distances of neighbors within distance maxr
    neigh_list=struc.get_all_neighbors(r=maxr,include_index=True)
    # for filtering distances from only species 1
    speclist=struc.indices_from_symbol(spec1)
    
    for icat in range(natoms): 
        # check for species
        if (icat in speclist):
            #print(neigh_list[icat])
            n_num = len(neigh_list[icat]) # num neighbors of atom icat
            
            n0index = icat
            neighind =[]  #indices of neighbors
            for neigh in range(n_num):
                neighind.append(neigh_list[icat][neigh][2])
                
            #find all pairs of neighbors by indices using itertools.combinations
            pairs=combinations(neighind,2)
            # find angle with previously defined central atom
            # and all possible pairs of nearest neighbor atoms
            for p in pairs:
                n1index = p[0]
                n2index = p[1]
                angle = getangle(struc,n1index,n0index,n2index)
                
                #for some reason the order of index  matters...., kluge fix
                tmpangle = getangle(struc,n2index,n0index,n1index)
                if tmpangle > angle:
                    angle = tmpangle
                    
                if (angle <= maxangle) and (angle>= minangle):
                    ang.append(angle);
        
    return ang 

#====================================================================================#

#------------------- Species specific distribution @ surface ------------------------#
def surf_spec_bond_len(struc,maxr,spec1,spec2,n=2):
    """Finding number of nearest Neighbors of distance r, version 2, with species distinction
    search for atoms within r radius (I think in units of lattice vectors) of each atomic site
    and bin the distances in histogram
  
    adapted for surface, look for atoms within n layers of the surface 

    get unique number of bond lengths
    give maxr slightly more than 2nd nn distance, 
    identify by element species 
 
    Input
     struc = structure read from pymatgen
     natoms = number of atoms, integer
     maxr = radius for searching distances for each site i, float;units of Ang.
     spec1, spec2 = identity of species associated with distance, chemical symbol; string
     n = integer, number of layers from surface to consider; +/- z
    Output
     rad: list of bond distances between spec1 and spec2
     surfatoms: list of Atoms included in histogram
    """
    
    rad = []
    natoms = struc.num_sites

    # get distances of neighbors within distance maxr
    neigh_list=struc.get_all_neighbors(r=maxr,include_index=False)
    # for filtering distances from only species 1
    speclist=np.array(struc.indices_from_symbol(spec1))

    #assumes slab constructed long-way along z
    #find z-cutoff for atoms of species 1 considered at surface 
    zcoords = np.zeros(len(speclist))
    for ispec in np.arange(len(speclist)):
      spec = speclist[ispec]
      zcoords[ispec] = struc.sites[spec].coords[2]

    # +/- z, average of n and n+1 inner layer z coordinates 
    uzcoords = np.unique(zcoords.round(decimals=2))
    uzcoords.sort()
    mz = np.average([uzcoords[n],uzcoords[n-1]])
    pz = np.average([uzcoords[-n],uzcoords[-n-1]])

    surfatoms = []
    for icat in range(natoms): 
      # check within +/- z of surface  
      z = struc.sites[icat].coords[2]
      if (z > pz) or (z < mz):
          surfatoms.append(struc.sites[icat])
          # check for species
          if (icat in speclist):
              #print len(neigh_list[itr])
              n_num = len(neigh_list[icat])
              for ineigh in range(n_num):
                  point=neigh_list[icat][ineigh]
                  #print point
                  #check if finding distance to spec 2 
                  tmpspec2=point[0].as_dict()['species'][0]['element']
                  if (tmpspec2 == spec2):
                      dist = point[1]
                      rad.append(dist);
                      #print spec1,spec2,tmpspec2,dist
    return rad, surfatoms
 
def surf_spec_bond_angle(struc,maxr,spec1,spec2,minangle=0.0,maxangle=180.0,n=2):
    """Finding bond angle distribution with spec 1 as the center atom, crystalline with definite 
        coordination number
    search for atoms within r radius (I think in units of lattice vectors) of each atomic site
    and bin the distances in histogram
 
    adapted for looking at atoms within n layers from surface 

    Inputs
     struc = structure read from pymatgen
     maxr = radius for searching distances for each site i, float; I think in units of Ang.
     spec1, spec2 = identity of species associated with distance; string
     minangle, maxangle = bounds of angles looking for
     n = integer, number of layers on both sides of slab to include
    Outputs
     ang = list of angles between for bonds of spec2-spec1-spec2
    """
    
    #from . import getangle # needed? 
  
    ang = []
    natoms = struc.num_sites
    
    # get distances of neighbors within distance maxr
    neigh_list=struc.get_all_neighbors(r=maxr,include_index=True)
    # for filtering distances from only species 1
    speclist=struc.indices_from_symbol(spec1)

    #assumes slab constructed long-way along z
    #find z-cutoff for atoms of species 1 considered at surface 
    zcoords = np.zeros(len(speclist))
    for ispec in np.arange(len(speclist)):
      spec = speclist[ispec]
      zcoords[ispec] = struc.sites[spec].coords[2]

    # +z, average of n and n+1 inner layer z coordinates 
    uzcoords = np.unique(zcoords.round(decimals=2))
    uzcoords.sort()
    mz = np.average([uzcoords[n],uzcoords[n-1]])
    pz = np.average([uzcoords[-n],uzcoords[-n-1]])

    
    for icat in range(natoms): 
       # check within +/- z of surface  
       z = struc.sites[icat].coords[2]
       if (z > pz) or (z < mz):
           # check for species
           if (icat in speclist):
               #print(neigh_list[icat])
               n_num = len(neigh_list[icat]) # num neighbors of atom icat
               
               n0index = icat
               neighind =[]  #indices of neighbors
               for neigh in range(n_num):
                   neighind.append(neigh_list[icat][neigh][2])
                   
               #find all pairs of neighbors by indices using itertools.combinations
               pairs=combinations(neighind,2)
               # find angle with previously defined central atom
               # and all possible pairs of nearest neighbor atoms
               for p in pairs:
                   n1index = p[0]
                   n2index = p[1]
                   angle = getangle(struc,n1index,n0index,n2index)
                   
                   #for some reason the order of index  matters...., kluge fix
                   tmpangle = getangle(struc,n2index,n0index,n1index)
                   if tmpangle > angle:
                       angle = tmpangle
                       
                   if (angle <= maxangle) and (angle>= minangle):
                       ang.append(angle);
        
    return ang 
