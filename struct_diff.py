#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: wwwennie

 returns metrics for characterizing difference between two given structures
 e.g., between relaxation runs

 in order for results to be meanginful, ensure that the index of the listed
 atoms is the same between structures

 based on some pymatgen functionality, so sensitive to atomic species
 but not occupation, spin
-----------------------------------------------------------
"""

import pymatgen as pmg
import numpy as np

# maybe useful later: 
#   pymatgen.analysis.structure_matcher; class StructureMatcher
#   http://pymatgen.org/pymatgen.analysis.structure_matcher.html


#==========================================================================

#--------------------- Compare relaxation -------------------------------#

def coord_diff(struc1,struc2,scale=0.5):
  """ Calculate difference between structures, in terms of atomic coordinates
      In calculating a translation vector
      assumes the the index of atoms remains same between compared structures

      Inputs:
        struc1, struc2 = pymatgen structures; everything is referenced to struc1, 
                         the 'original' structure
        scale = how much to scale the magnitude of change (for visualization purposes)
 
      Outputs: 
         uvw = list of vectors describing change in atomic positions
  """

  coord1 = struc1.cart_coords
  coord2 = struc2.cart_coords

  uvw = np.subtract(coord2,coord1)

  return coord1,uvw*scale 
  
  
