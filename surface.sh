#!/usr/bin/env bash

"""
--------------------------------------------------------------
 created: wwwennie
 Wrapper for making surface structure from QE -> pymatgen

 Usage:    surface.sh <prefix>.in <miller index> slabsize vacsize
    e.g.,  surface.sh <prefix>.in '1,1,1' 5 5 
             for all terminations of specific miller index
    e.g.,  surface.sh <prefix>.in 1 5 5 
             for all terminations of all indices with 
             max index of miller index
  where slabsize and vacsize is minimum slab and vacuum size in Ang
  see pymatgen documentation for more info

 assumes Xcrysden is installed and pwi2xsf & pwo2xsf are in PATH
 depends on qe2pmg.sh
--------------------------------------------------------------
"""

ubin="/home/wwwennie/wwwennie@uchicago.edu/bin/struct/"
filname=$1
filout=$(echo "${filname%%.*}")


# check input or output
if [[ $filname == *.in ]]; then
  pwi2xsf $filname > tmp.xsf 
else
  pwo2xsf $filname > tmp.xsf 
fi 

# make compatible with pymatgen readin
${ubin}/qe2pmg.sh tmp.xsf 
mv out-tmp.xsf $filout-pmg.xsf

# make surface of specified miller index
echo "Creating surfaces"
python ${ubin}/surface.py $filout-pmg.xsf $2 $3 $4

# convert back to QE input
echo "==== *** Reminder: Convert from cif to QE input *** ===="
#${ubin}/cif2qe.sh super-$size-$filout-pmg.cif > super-$size-$filout.in

rm $filout-pmg.xsf
#rm super-$size-$filout-pmg.cif
