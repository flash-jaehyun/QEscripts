#!/usr/bin/env bash

# created: wwwennie
# wrapper for converting QE -> XSF -> pymatgen
# pymatgen (2018.8.7)  is unable to handle element symbols, only atomic numbers
# this is the wrapper for mapping to element symbols to atomic numbers

# usage: ./qe2pmg.sh <file name>.xsf

# thought about hash tables, but not sed and awk not compatible
# hackneyed way to replace all elements
# up to Bi, skipped lanthanides

declare -A periodictable=( ["H"]="1" ["He"]="2" ["Li"]="3" ["Be"]="4" ["B"]="5" \
                        ["C"]="6" ["N"]="7" ["O"]="8" ["F"]="9" ["Ne"]="10"\
			["Na"]="11" ["Mg"]="12" ["Al"]="13" ["Si"]="14" ["P"]="15" \
			["S"]="16" ["Cl"]="17" ["Ar"]="18" ["K"]="19" ["Ca"]="20" \
			["Sc"]="21" ["Ti"]="22" ["V"]="23" ["Cr"]="24" ["Mn"]="25" \
			["Fe"]="26" ["Co"]="27" ["Ni"]="28" ["Cu"]="29" ["Zn"]="30" \
			["Ga"]="31" ["Ge"]="32" ["As"]="33" ["Se"]="34" ["Br"]="35" \
			["Kr"]="36" ["Rb"]="37" ["Sr"]="38" ["Y"]="39" ["Zr"]="40" \
			["Ni"]="41" ["Mo"]="42" ["Tc"]="43" ["Ru"]="44" ["Rh"]="45" \
			["Pd"]="46" ["Ag"]="47" ["Cd"]="48" ["In"]="49" ["Sn"]="50" \
			["Sb"]="51" ["Te"]="52" ["I"]="53" ["Xe"]="54" ["Cs"]="55" \
			["Ba"]="56" ["La"]="57" \
			["Hf"]="72" ["Ta"]="73" ["W"]="74" ["Re"]="75"\
			["Os"]="76" ["Ir"]="77" ["Pt"]="78" ["Au"]="79" ["Hg"]="80"\
			["Tl"]="81" ["Pb"]="82" ["Bi"]="83")
#echo "${periodictable[Bi]}"

infil=$1

# output head of xsf, it is unchanged
sed '/PRIMCOORD/q' $infil > tmp0

# isolate atomic coordinates
sed -n -e '/PRIMCOORD/,$p' $infil > tmp

# extract element names, find unique set
elements=($(awk 'NR>2 {print $1}' tmp))
natoms=${#elements[@]}
uniq=($(awk 'NR>2 {print $1}' tmp | sort -u))
nuniq=${#uniq[@]}

# replace element name with symbol
for el in "${uniq[@]}"
do
	:
	sed -i "s/$el/${periodictable[$el]}/g" tmp
done

# account for needed headers
head=$((natoms+1))
tail -n $head tmp >> tmp2

# tack onto output file
cat tmp0 tmp2 >> out-$infil

rm tmp*
rm pwi2xsf.xsf_out
rm pw.*
