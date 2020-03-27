#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: wwwennie
 
 Quantum epsresso input to cif output
---------------------------------------------------------
"""

from ase.io import read,write
from ase.build import molecule, add_adsorbate
from pymatgen.io.cif import CifWriter

# from struct-<>.py
from struct_import import *
from struct_defects import *
from struct_orient import orient

# import structure into ASE
filname = "relax.in"
bivo4 = struct_import(filname, to_pymatgen=False) 
write("relax.cif",bivo4,format="cif")

