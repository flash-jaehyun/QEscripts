#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: wwwennie
 
 mini demo for using series of struct-<>.py files for 
  analysing structures using pymatgen

 TODO: move to Jupyter notebook for plotting 
  
 here, I'm just testing functionality as it gets scripted

---------------------------------------------------------
"""

import subprocess
from pymatgen import Lattice, Structure, Site
from pymatgen.io.xcrysden import XSF

# from struct-<>.py
from struct_import import struct_import

# import structure into pymatgen
filname="tet.in"
#ubin="/home/wwwennie/wwwennie@uchicago.edu/bin/struct/"
strucbin="/home/wwwennie/wwwennie@uchicago.edu/bin/pos_tmp/"

filname = strucbin+filname
bivo4 = struct_import(filname)

##### useful printing statements #####
#(dist,transvec) = bivo4.sites[16].distance_and_image(bivo4.sites[4])
#print(len(bivo4.sites))
#print(bivo4.sites)
#print(bivo4.frac_coords)
#print(bivo4.sites[1])
#print(bivo4.sites[1].coords)
#print(bivo4.lattice.abc)
#print(bivo4.sites[1].frac_coords)
#print(bivo4.lattice.matrix)

