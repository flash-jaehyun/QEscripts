#!/usr/bin/env bash

"""
---------------------------------------------------------------
 created: wwwennie@gmail.com
 wrapper for making supercell structure from QE -> pymatgen
 pymatgen (2018.8.7)  is unable to handle element symbols, only atomic numbers
                      and is not quite working with QE input/output files

 usage: supercell.sh <prefix>.in '2,3,4'
        supercell.sh <prefix>.out 2
        supercell.sh <prefix>.out '2,1,0' '0,3,0' '0,0,1'

 depends on supercell.py
 see pymatgen documentation for more information
 
 assumes Xcrysden is installed and pwi2xsf & pwo2xsf are in PATH
 depends on qe2pmg.sh
---------------------------------------------------------------
"""

ubin="/home/wwwennie/wwwennie@uchicago.edu/bin/QEscripts/" #TODO: change this
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

# make supercell
# variable input arguments allowed
if [[ $# -eq 2 ]]; then
  python $ubin/supercell.py $filout-pmg.xsf $2
 
  if [[ ${#2} -eq 1 ]]; then 
    size=$2$2$2
  else
    size=$(echo "$2" | sed -r 's/[',']+//g')
  fi
 
elif [[ $# -eq 4 ]]; then
  python $ubin/supercell.py $filout-pmg.xsf $2 $3 $4

  l=$(echo "$2" | sed -r 's/[',']+//g')
  m=$(echo "$3" | sed -r 's/[',']+//g')
  n=$(echo "$4" | sed -r 's/[',']+//g')
  
  size=$l'x'$m'x'$n
fi 


# convert back to QE input
echo "Converting from cif to QE input"
${ubin}/cif2qe.sh super-$size-$filout-pmg.cif > super-$size-$filout.in

rm $filout-pmg.xsf
#rm super-$size-$filout-pmg.cif
