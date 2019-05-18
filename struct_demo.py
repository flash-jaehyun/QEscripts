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
from pymatgen.io.cif import CifWriter

# from struct-<>.py
from struct_import import struct_import
from struct_bonds import getangle,cryst2cart
import struct_bonds as bonds
import struct_polyhedra as poly
import struct_plot as plot
import struct_diff 
from struct_defects import *

# import structure into pymatgen
strucbin="/home/wwwennie/wwwennie@uchicago.edu/bin/pos_tmp/"

#filname = strucbin+"tet.in"
#bivo4 = struct_import(ubin,filname)
#
#filname = strucbin+"pseudotet.in"  # volume expansion, not quite hydrostatic
#bivo4comp = struct_import(filname)

filname = strucbin+"Ov-s15v20.in"
bivo4surf = struct_import(filname)

########### polaron bias #########
struct = bias_polaron(bivo4surf,121,1.9,bias=1.05)
w = CifWriter(struct)
w.write_file("bias.cif")

########### surface structure histogram #######
#title="Bond length distribution"
#legend=["Bi-O","V-O"]
#len_BiO,surfBiatoms = bonds.surf_spec_bond_len(bivo4surf,2.8,'Bi','O',n=2)
#len_VO,surfVatoms = bonds.surf_spec_bond_len(bivo4surf,2.8,'V','O',n=2)
#plot.hist_plot([len_BiO,len_VO], title,legend,xlim=[1.6,2.65])
#
#title="Bond angle distribtion"
#legend=["O-Bi-O","O-V-O"]
#ang_BiO = bonds.surf_spec_bond_angle(bivo4,2.8,'Bi','O',n=2)
#ang_VO = bonds.surf_spec_bond_angle(bivo4,2.8,'V','O',n=2)
#plot.hist_plot([ang_BiO,ang_VO], title,legend,xlim=[65,165])
#
###### effect of relaxation #######
#coords,uvw = struct_diff.coord_diff(bivo4,bivo4comp,scale=0.5)
#print(coords)
#plot.quiver3d(coords,uvw)

###### plotting #########
#title="Bond angle distribtion"
#legend=["O-Bi-O","O-V-O"]
#ang_BiO = bonds.spec_bond_angle(bivo4,2.8,'Bi','O')
#ang_VO = bonds.spec_bond_angle(bivo4,2.8,'V','O')
#plot.hist_plot([ang_BiO,ang_VO], title,legend)
#
#title="Bond length distribution"
#legend=["Bi-O","V-O"]
#len_BiO = bonds.spec_bond_len(bivo4,2.8,'Bi','O')
#len_VO = bonds.spec_bond_len(bivo4,2.8,'V','O')
#plot.hist_plot([len_BiO,len_VO], title,legend)

###### polyhedron statistics ######
#angle_var = poly.bond_ang_var(bivo4,2.8,6)
#print(angle_var)
#quad_elong = poly.bond_quad_elong(bivo4,2.8,6)
#print(quad_elong)

#### coordination number #####
#cn = poly.get_all_cn(bivo4,2.7,'Bi')
#print(cn)

###### testing spec_bond_len #####
#bi_o = bonds.spec_bond_len(bivo4,4.8,'Bi','O')
#print(bi_o)
#
###### testing spec_bond_len #####
#o_bi_o = bonds.spec_bond_angle(bivo4,2.8,'Bi','O')
#print(o_bi_o)

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

###### structure analysis: testing getangl #####
# indexing based on python indexing
## tetragonal cell
#print("same cell", getangle(bivo4,34,5,53)) # same cell
#print("---")
#print("+x", getangle(bivo4,51,3,35)) #  +x cell, o-bi-o, 137.346 
#print("---")
#print("+y",getangle(bivo4,6,63,12)) # +y cell, v-o-bi, 122.121
#print("---")
#print("+z",getangle(bivo4,10,46,7)) # +z cell, bi-o-bi, 104.233
#print("---")
#print("-z",getangle(bivo4,7,46,10)) # -z cell

## monoclinic cell
#print("same cell", getangle(bivo4,40,17,38)) # same cell, 135.309
#print("---")
#print("+x", getangle(bivo4,68,18,69)) #  +x cell, o-bi-o, 135.309 
#print("---")
#print("+y",getangle(bivo4,10,49,1)) # +y cell, v-o-v, 107.325
#print("---")
#print("+z",getangle(bivo4,27,3,24)) # +z cell, o-v-o, 115.457
#print(getangle(bivo4,19,57,16)) # 104.310
#print("--small angles--")
#print(getangle(bivo4,39,24,3)) # 37.377, same cell
#print(getangle(bivo4,30,0,33)) # 40.842, +z
