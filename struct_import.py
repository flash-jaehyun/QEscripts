#!/usr/bin/env python3
"""
-----------------------------------------------------------
 created: wwwennie

  import QE structure from input or output file
   and import into Structure class for pymatgen
   manipulation in python
 
  see pymatgen documentation for more information
 
  There are several options
   1) through ASE interface, struct_import
   2) through XCRYSDEN scripts, struct_import_xcrys

  Inputs 
     ubin:   string, directory where qe2pmg.sh is 
     filname: QE file with atomic coordinates, assumes *.in or *.out
  Outputs
     structure:  pymatgen.structure
----------------------------------------------
"""

import subprocess
import os
from pymatgen import Lattice, Structure


def struct_import(filname,to_pymatgen=True):
   """ import structure from QE input/output into Pymatgen structure
       through ASE interface
  
       goal: make this module standalone independent of XCRYSDEN or anything outside of python
       intermediate: go through ASE
       eventually: adapt rest of struct module with ASE
  
       filname: input/output QE file containing structure
       ubin: (string) dummy variable
       returns pymatgen structure if to_pymatgen, otherwise ASE Atoms object
   """
   from ase.io import read
   
   print("# Reading atomic coordinates from: ", filname)
   # compatbility with QE input/output
   if filname.endswith(".in"):
      structure = read(filname,index=None,format="espresso-in")
   elif filname.endswith(".out"):
      structure = read(filname,index=-1,format="espresso-out")
   
   if to_pymatgen:
      from pymatgen.io.ase import AseAtomsAdaptor
      structure = AseAtomsAdaptor.get_structure(structure)

   return structure


def struct2ase(struct):
   """ convert pymatgen structure to ASE structure """
   from pymatgen.io.ase import AseAtomsAdaptor
   structure = AseAtomsAdaptor.get_structure(struct)
   return structure

def struct2pmg(struct):
   """ convert ASE structure to pymatgen structure 
       this is not natively implemented in ASE
   """
   from ase.io import write
   from pymatgen.io.cif import CifParser # this didn't give supercell structure
   from pymatgen.io.vasp import Poscar
   write("POSCAR.tmp",struct,format="vasp")
  
   # fix formatting of the outputted vasp format, which has misplaced species labels
   # this is a hacky way to insert it in the right line
   f = open("POSCAR.tmp","r")
   lines = []
   for i,line in enumerate(f):
      if i == 0:
        species_line = line
      elif i == 5:
        lines.append(species_line)
      lines.append(line)
   f.close()
   with open("POSCAR",'w') as f:
     f.writelines(lines)  


   poscar = Poscar.from_file("POSCAR")
   structure = poscar.structure
   subprocess.call(["rm", "POSCAR","POSCAR.tmp"])

   return structure   

def struct_import_xcrys(filname,ubin=""):
   """ import structure from QE input/output into Pymatgen structure
       essentially creates intermediary file *.xsf that both programs can parse
  
       ubin: directory containing qe2pmg.sh
       filname: input/output QE file containing structure

       unlike struct_import, can also handle PWSCF output file
       but requires installation of XCRYSDEN
 
       NOTE: not as portable
       TODO: make standalone; 
         roadblock: pw*2xsf executuable in XCRYSDEN_LIB_BINDIR is binary, so cannot tell
                    what other dependencies there are
       
   """
   #formatting outputname
   outname=filname.split('.')
   outname=outname[0].split('/')
   outname=outname[-1:]

   print("NOTE: XCRYSDEN INSTALLATION REQUIRED")
   print("# Reading atomic coordinates from: ", filname)
   # compatbility with QE input/output
   if filname.endswith(".in"):
       with open("out.xsf","w") as outfile:
           subprocess.check_call(["pwi2xsf",filname],stdout=outfile)
   elif filname.endswith(".out"):
       with open("out.xsf","w") as outfile:
           subprocess.check_call(["pwo2xsf",filname],stdout=outfile)

   # no check_call, since it returns a non-zero, which raises error...
   subprocess.call([str(ubin)+"qe2pmg.sh", "out.xsf"])
   subprocess.check_call(["mv","out-out.xsf", "out-pmg.xsf"])
   
   # Read structure from PWscf file, xsf format
   structure = Structure.from_file("out-pmg.xsf")
   
   subprocess.call(["rm", "out.xsf","out-pmg.xsf"])
   subprocess.call(["rm","-r","__pycache__"])

   return structure

def struct_import_pymatgen(filename,ubin=""):
   """ import structure from QE input/output into Pymatgen structure
       adapted from pymatgen native PWio method
             pymatgen.io.pwscf from_file and from_string methods
  
       DOES NOT WORK ATM
  
       filname: (string) input/output QE file containing structure
       ubin: (string) dummy variable; for consistency with remaining scripts
   """
   import re
   from monty.io import zopen
   from monty.re import regrep 
   from pymatgen.util.io_utils import clean_lines
   from pymatgen.io.pwscf import PWInput

   print("# Reading atomic coordinates from: ", filename)
   with zopen(filename, "rt") as f:
            string = f.read()

   lines = list(clean_lines(string.splitlines()))

   print(lines)
   def input_mode(line):
       if line[0] == "&":
           return ("sections", line[1:].lower())
       elif "ATOMIC_SPECIES" in line:
           return ("pseudo", )
       elif "K_POINTS" in line:
           return ("kpoints", line.split("{")[1][:-1])
       elif "CELL_PARAMETERS" in line or "ATOMIC_POSITIONS" in line:
           return ("structure", line.split("{")[1][:-1])
       elif line == "/":
           return None
       else:
           return mode

   sections = {"control": {}, "system": {}, "electrons": {}, 
               "ions": {}, "cell":{}}
   pseudo = {}
   pseudo_index = 0
   lattice = []
   species = []
   coords = []
   structure = None
   site_properties = {"pseudo":[]}
   mode = None
   for line in lines:
       mode = input_mode(line)
       if mode == None:
           pass
       elif mode[0] == "sections":
           section = mode[1]
           m = re.match(r'(\w+)\(?(\d*?)\)?\s*=\s*(.*)', line)
           if m:
               key = m.group(1).strip()
               key_ = m.group(2).strip()
               val = m.group(3).strip()
               if key_ != "":
                   if sections[section].get(key, None) == None:
                       val_ = [0.0]*20 # MAX NTYP DEFINITION
                       val_[int(key_)-1] = PWInput.proc_val(key, val)
                       sections[section][key] = val_

                       site_properties[key] = []
                   else:
                       sections[section][key][int(key_)-1] = PWInput.proc_val(key, val) 
               else:
                   sections[section][key] = PWInput.proc_val(key, val)

#       elif mode[0] == "pseudo":
#           m = re.match(r'(\w+)\s+(\d*.\d*)\s+(.*)', line)
#           if m:
#               pseudo[m.group(1).strip()] = {}
#               pseudo[m.group(1).strip()]["index"] = pseudo_index
#               pseudo[m.group(1).strip()]["pseudopot"] = m.group(3).strip()
#               pseudo_index += 1
#       elif mode[0] == "kpoints":
#           m = re.match(r'(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', line)
#           if m:
#               kpoints_grid = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
#               kpoints_shift = (int(m.group(4)), int(m.group(5)), int(m.group(6)))
#           else:
#               kpoints_mode = mode[1]
       elif mode[0] == "structure":
           m_l = re.match(r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)', line)
           m_p = re.match(r'(\w+)\s+(-?\d+\.\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)', line)
           if m_l:
               lattice += [ float(m_l.group(1)), float(m_l.group(2)), float(m_l.group(3)) ]
           elif m_p:
               site_properties["pseudo"].append(pseudo[m_p.group(1)]["pseudopot"])
               species += [pseudo[m_p.group(1)]["pseudopot"].split(".")[0]]
               coords += [[float(m_p.group(2)), float(m_p.group(3)), float(m_p.group(4))]]

               for k, v in site_properties.items():
                   if k != "pseudo":
                       site_properties[k].append(sections['system'][k][pseudo[m_p.group(1)]["index"]])
           if mode[1] == "angstrom":
               coords_are_cartesian = True
           elif mode[1] == "crystal":
               coords_are_cartesian = False

   structure = Structure(Lattice(lattice), species, coords, 
                              coords_are_cartesian=coords_are_cartesian,
                              site_properties=site_properties)
   return structure


