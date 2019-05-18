#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: wwwennie

  find distributions of bond lengths, bond angles of polyhedra 


  see pymatgen documentation for more information

----------------------------------------------
"""

import numpy as np
from itertools import combinations
from struct_bonds import getangle_sites

#------------------------------- Helper functions --------------------------------#
def get_polyhedron(struc,maxr,siteindex):
    """ Finding a particular polyhedra surrounding spec1 for some distance maxr
        set maxr to just slightly above nn distance

    Input
      struc = structure read from pymatgen
      maxr = radius of search; in ang.
      siteindex = string; chemical name of species of interest
   
    Output 
      poly = list of sites corresponding to polyhedra
             of form [(site,dist)...]
      cn = coordination number of polyhedra
    """

    site = struc.sites[siteindex]
    poly = struc.get_neighbors(site,maxr)
    cn = len(poly)

    if (cn <= 2):
      print("You didn't choose a polyhedron! CN <= 2")
      return None

    return poly, cn


def get_all_cn(struc,maxr,spec1):
    """Finding coordination number of all atoms of species 1 within distance maxr

     Input:
      struc = structure read from pymatgen
      maxr = radius for searching distances for each site i, float; I think in units of Ang.
      spec1, spec2 = identity of species associated with distance; string

     Ouput:
      cn = list of coordination numbers of spec1 within radius r
    """
    cn = []
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
            cn.append(n_num);
                    #print spec1,spec2,tmpspec2,dist

 
    return cn

def bond_quad_elong(struc,maxr,siteindex,ideal_len=0.0):
    """ Find quadratic elongation for polyhedra surrounding siteindex
        inspired by: 10.1126/science.172.3983.567, specific to octehedra, but 
           probably ok for  

    Input
      struc = structure read from pymatgen
      maxr = float; radius of search; in ang.
      siteindex = string; chemical name of species of interest
      ideal_len = float; angstroms; ideal bond length; if not given, defaults 
                           to avg of bond lengths 
   
    Output 
     qelong = float; quadratic elongation
    """
 
    poly,cn = get_polyhedron(struc,maxr,siteindex)
    poly_sites = [i[0] for i in poly]
    dist = [i[1] for i in poly]

    if (ideal_len == 0.0):
       ideal_len = np.mean(dist)
   
    qelong = np.mean(np.divide(dist,ideal_len))

    return qelong
  

def bond_ang_var(struc,maxr,siteindex):
    """ Find bond angle variance for polyhedra surrounding siteindex
        inspired by: 10.1126/science.172.3983.567, which is specific to octahedra
        here, the variance is unbiased 

    Input
      struc = structure read from pymatgen
      maxr = radius of search; in ang.
      siteindex = string; chemical name of species of interest
   
    Output 
     var = float; bond angle variance

    TODO: 
      - way to filter out bonds that are not between adjacent corners....
    """
 
    poly,cn = get_polyhedron(struc,maxr,siteindex)
    poly_sites = [i[0] for i in poly]


    ang = []
    # find all possible combinations of nearest neighbors 
    pairs = combinations(poly_sites,2)
    site2 = struc.sites[siteindex]
    for p in pairs:
      site1 = p[0]
      site3 = p[1]
      angle = getangle_sites(struc,site1,site2,site3)
      ang.append(angle)

    # for limited set of polyhedra: tetrahedra, octahedra, optohedra (?)
    #   ideal, undistorted angle between corner-center atom-corner
    if (cn == 4): 
        ideal = 109.5*np.ones(len(ang))
    elif (cn == 6):
        ideal = 90.0*np.ones(len(ang))
    elif (cn == 8):
        ideal = 76.0*np.ones(len(ang)) # based on average angle in XeFe8-
    else:
        ideal = np.mean(ang)
    
    #print(struc.sites[siteindex],ideal,ang, np.mean(ang))
    #var =  np.sum(abs(ang - np.mean(ang))**2)/(float(len(ang))-1.0)
    var =  np.sum(abs(ang - ideal)**2)/(float(len(ang))-1.0)
    print("Take with grain of salt; this includes ALL possible angles, included non-adjacent neighbors")
    return var
  
