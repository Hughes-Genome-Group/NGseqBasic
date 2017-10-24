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

setPloidyPath(){
    
ploidyPath="UNDETERMINED"
    
for g in $( seq 0 $((${#genomesWhichHaveBlacklist[@]}-1)) ); do

# echo ${supportedGenomes[$g]}

if [ "${genomesWhichHaveBlacklist[$g]}" == "${genomeName}" ]; then
    ploidyPath="${BLACKLIST[$g]}"
fi

done 
    
if [ "${ploidyPath}" == "UNDETERMINED" ]; then
  echo
  echo "NOTE !! Genome build " ${genomeName} " is not supported in BLACKLIST FILTERING - turning blacklist filtering off !"
  echo
  echo  >&2
  echo "NOTE !! Genome build " ${genomeName} " is not supported in BLACKLIST FILTERING - turning blacklist filtering off !"  >&2
  echo  >&2
  
  weHavePloidyFile=0;

else 
  weHavePloidyFile=1;
fi

# Check that the file exists, and informing user that the path was set ..
if [ "${weHavePloidyFile}" -eq 1 ] ; then

if [ ! -e "${ploidyPath}" ] || [ ! -r "${ploidyPath}" ] || [ ! -s "${ploidyPath}" ]; then
  echo "Blacklisted regions for genome ${genomeName} file ${ploidyPath} not found or empty file - aborting !"  >&2
  exit 1     
fi

echo
echo "Set BLACKLIST FILTERING file : ${ploidyPath}"

fi

}

