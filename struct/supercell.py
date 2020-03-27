#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: wwwennie
 Wrapper for making supercell from QE -> pymatgen

 Use case: python supercell.py <QE input/output> int int int
  e.g., python supercell.py bivo4.xsf 2,3,4
      for 2x3x4 supercell 
  e.g., python supercell.py bivo4.xsf 2
      for 2x2x2 supercell 
  e.g., python supercell bivo4.xsf 2,1,0 0,3,0 0,0,1
      for new structure with a' = 2a+b, b' = 3b, c'= c

  see pymatgen documentation for more information

  assumes Xcrysden is installed and pwi2xsf & pwo2xsf are in PATH
  depends on qe2pmg.sh
----------------------------------------------
"""

import sys
from pymatgen import Lattice, Structure
from pymatgen.io.xcrysden import XSF
from pymatgen.io.cif import CifWriter

# Input parameters, following behavior of make_supercell in pymatgen
filname=sys.argv[1]
outname=filname.split('.')
outname=outname[0].split('/')
outname=outname[-1:]

print("Reading atomic coordinates from: ", filname)

# formatting string for output file name
if len(sys.argv) == 3:
    sp=sys.argv[2].split(',')
    if len(sys.argv[2]) == 1:
      size=sys.argv[2]*3
    else:
      size=sys.argv[2].replace(',','')
elif len(sys.argv) == 5:
    a=sys.argv[2].split(',')
    b=sys.argv[3].split(',')
    c=sys.argv[4].split(',')
    sp =[a,b,c]
 
    l=sys.argv[2].replace(',','')
    m=sys.argv[3].replace(',','')
    n=sys.argv[4].replace(',','')
    size=l+"x"+m+"x"+n
else:
    print("Check input formatting")

## compatbility with QE input/output
#if filname.endswith(".in"):
#    with open("fil.xsf","w") as outfile:
#        subprocess.check_call(["pwi2xsf",filname],stdout=outfile)
#    filname="fil.xsf"
#elif filname.endswith(".out"):
#    subprocess.check_call(["pwo2xsf",filname," > fil.xsf"])
#    filname="fil.xsf"

# Read structure from PWscf file, xsf format
structure = Structure.from_file(filname)

# make supercell
structure.make_supercell(sp,to_unit_cell=True)

# output to cif 
w = CifWriter(structure)
w.write_file("super-"+size+"-"+outname[0]+".cif")
