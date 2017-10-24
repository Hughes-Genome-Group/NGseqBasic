#!/bin/bash

##########################################################################
# Copyright 2016, Jelena Telenius (jelena.telenius@imm.ox.ac.uk)         #
#                                                                        #
# This file is part of NGseqBasic .                                      #
#                                                                        #
# NGseqBasic is free software: you can redistribute it and/or modify     #
# it under the terms of the GNU General Public License as published by   #
# the Free Software Foundation, either version 3 of the License, or      #
# (at your option) any later version.                                    #
#                                                                        #
# NGseqBasic is distributed in the hope that it will be useful,          #
# but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
# GNU General Public License for more details.                           #
#                                                                        #
# You should have received a copy of the GNU General Public License      #
# along with NGseqBasic.  If not, see <http://www.gnu.org/licenses/>.    #
##########################################################################

setBOWTIEgenomeSizes(){
    
bowtieGenomeBuild="UNDETERMINED"

#-----------Genome-sizes-for-bowtie-commands----------------------------------------------

if [ "${BOWTIE}" -eq 1 ]
then
    
    
for g in $( seq 0 $((${#supportedGenomes[@]}-1)) ); do
    
# echo ${supportedGenomes[$g]}

if [ "${supportedGenomes[$g]}" == "${genomeName}" ]; then
    bowtieGenomeBuild="${BOWTIE1[$g]}"
fi

done  

fi

#------------------------------------------------

if [ "${BOWTIE}" -eq 2 ]
then
    
for g in $( seq 0 $((${#supportedGenomes[@]}-1)) ); do
    
# echo ${supportedGenomes[$g]}

if [ "${supportedGenomes[$g]}" == "${genomeName}" ]; then
    bowtieGenomeBuild="${BOWTIE2[$g]}"
fi

done  

fi  

#------------------------------------------------

# Check that it got set ..

if [ "${bowtieGenomeBuild}" == "UNDETERMINED" ]; then 
  echo "Genome build " ${genomeName} " is not supported - aborting !"  >&2
  exit 1 
fi    

# Check that at least one index file exists..

TEMPcount=$(($( ls -1 ${bowtieGenomeBuild}* | grep -c "" )))

if [ "${TEMPcount}" -eq 0 ]; then

  echo "Bowtie genome build for ${genomeName} : no index files ${bowtieGenomeBuild} found - aborting !"  >&2
  exit 1     
fi

echo
echo "Genome ${genomeName} .  Set BOWTIE index directory and basename : ${bowtieGenomeBuild}"

}



setUCSCgenomeSizes(){
    
genomeBuild="UNDETERMINED"
    
for g in $( seq 0 $((${#supportedGenomes[@]}-1)) ); do
    
# echo ${supportedGenomes[$g]}

if [ "${supportedGenomes[$g]}" == "${genomeName}" ]; then
    genomeBuild="${UCSC[$g]}"
fi

done 
    
if [ "${genomebuild}" == "UNDETERMINED" ]; then 
  echo "Genome build " ${genomeName} " is not supported - aborting !"  >&2
  exit 1 
fi

# Check that the file exists..
if [ ! -e "${genomeBuild}" ] || [ ! -r "${genomeBuild}" ] || [ ! -s "${genomeBuild}" ]; then
  echo "Genome build ${genomeName} file ${genomeBuild} not found or empty file - aborting !"  >&2
  exit 1     
fi

echo
echo "Genome ${genomeName} . Set UCSC genome sizes file : ${genomeBuild}"
echo

}



