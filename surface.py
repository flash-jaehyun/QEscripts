#!/usr/bin/env python3

"""
---------------------------------------------------------------
 created: wwwennie
 Wrapper for making surface structure from QE -> pymatgen

 Usage:    surface.py <prefix>.in <miller index> slabsize vacsize
    e.g.,  surface.py <prefix>.in '1,1,1' 5 5 
             for all terminations of specific miller index
    e.g.,  surface.py <prefix>.in 1 5 5 
             for all terminations of all indices with 
             max index of miller index
  where slabsize and vacsize is minimum slab and vacuum size in Ang
  see pymatgen documentation for more info

 assumes Xcrysden is installed and pwi2xsf & pwo2xsf are in PATH
 depends on qe2pmg.sh
------------------------------------------------------------------
"""

import sys
from pymatgen.core.surface import SlabGenerator,generate_all_slabs,Structure,Lattice
from pymatgen.io.xcrysden import XSF
from pymatgen.io.cif import CifWriter

# Input parameters, following behavior of SlabGenerator in pymatgen 
filname=sys.argv[1]
outname=filname.split('.')
outname=outname[0].split('/')
outname=outname[-1:]

print("Reading atomic coordinates from: ", filname)

# remaining input parameters
miller=sys.argv[2].split(',')
slabsize=int(sys.argv[3])
vacsize=int(sys.argv[4])

# Read structure from PWscf file, xsf format
structure = Structure.from_file(filname)

# make slab
if len(miller) == 3:
  miller=tuple(map(int,miller)) # no list in py2
  slabgen = SlabGenerator(structure,miller,slabsize,vacsize, center_slab=True)
  slabs = slabgen.get_slabs()
  print("Generated %d slabs in %s direction"%(len(slabs),str(miller),))

elif len(miller) == 1:
  miller=int(miller[0])
  slabs = generate_all_slabs(structure,miller,slabsize,vacsize,center_slab=True,symmetrize=True) 
  #slabs = generate_all_slabs(structure,miller,slabsize,vacsize,center_slab=True,symmetrize=False) 
  print("%s unique slab structures with max Miller index of %d"%(len(slabs),miller))

size="s"+str(slabsize)+"v"+str(vacsize)

counter = 1
for slab in slabs:
  print("Structure:", counter, slab.miller_index)
  print("Symmetric?", slab.is_symmetric())
  print("Polar?", slab.is_polar())
  counter = counter + 1

# output to cif 
counter=1
for slab in slabs:
  w = CifWriter(slab)
  index=''.join(map(str,slab.miller_index))
  w.write_file("slab-"+size+"-"+index+"-"+str(counter)+"."+outname[0]+".cif")
  counter = counter +1
