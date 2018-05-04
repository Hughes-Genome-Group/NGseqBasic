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

setPublicLocations(){

# #############################################################################

# This is the CONFIGURATION FILE to set up your SERVER ADDRESS and public data area
#  ( conf/serverAddressAndPublicDiskSetup.sh )

# Fill the name convention of your:

# - server type (ftp, http ..)
# - server address 
# - public disk space name convention
#    (if the disk address is not the same as the html address sub-path)

# As given in below examples

# ----------------------------------------------

# PUBLIC DATA - FOR UCSC VISUALISATION

# All visualisation data will be saved in a publicly available disk space

# If PIPE_hubbing.txt is used - ALL FILES are saved directly in the publicly available disk space.
# If PIPE_hubbingSymbolic.txt is used, symbolic links are generated for the bigwig files instead of storing the actual files.

# ----------------------------------------------

# The HTTP (or ftp) server address, and public data folder structure

SERVERTYPE="http"
SERVERADDRESS="userweb.molbiol.ox.ac.uk"

# These together set the server address to be :
# http://userweb.molbiol.ox.ac.uk

# Now, it assumes that disk area location
# /public/datafile.txt

# Is located in the server like this :
# http://userweb.molbiol.ox.ac.uk/public/datafile.txt

# If this is not the case - set up the OPTIONAL SETTINGS below
# (to finetune the relationship between disk space and server address)

# ----------------------------------------------

# OPTIONAL SETTINGS (if your disk area path differs from your server address sub-path) :

# You can delete or add stuff from the path, to match your data area design, with these parameters :

REMOVEfromPUBLICFILEPATH=""

# For example, your disk area location is :
# /upper/folder/public/datafile.txt

# And your publicly available file in server :
# http://userweb.molbiol.ox.ac.uk/public/datafile.txt

# To reach that, you would need to set :
# REMOVEfromPUBLICFILEPATH="/upper/folder"

ADDtoPUBLICFILEPATH=""

# For example, your disk area location is :
# /public/datafile.txt

# Publicly available file in server :
# http://userweb.molbiol.ox.ac.uk/home/folder/public/datafile.txt

# To reach that, you would need to set :
# ADDtoPUBLICFILEPATH="/home/folder"

tobeREPLACEDinPUBLICFILEPATH=""
REPLACEwithThisInPUBLICFILEPATH="" # This can be empty : then the 'tobeREPLACEDinPUBLICFILEPATH' will be replaced with '' (just removed)

# For example, your disk area location is :
# /public/sub/folder/datafile.txt

# And your publicly available file in server :
# http://userweb.molbiol.ox.ac.uk/public/another/folder/structure/datafile.txt

# To reach that, you would need to set :
# tobeREPLACEDinPUBLICFILEPATH="/sub/folder"
# REPLACEwithThisInPUBLICFILEPATH="/another/folder/structure"


# THE ORDER THE ABOVE TAKE EFFECT :

# 1) REMOVE from public path
# 2) ADD to public path
# 3) REPLACE in public path

# Note !! - the (3) is vulnerable to repeating of the same subfolder structure.
# It only will find and replace the LEFTmost of the occurrences of the file path it searches for.

# You can finetune the subroutine parsePublicLocations() which does this parsing :
# it is given below, in the end of this config file !

}

# -------------------------------------------
# The parser of server name and public file location.
# The parser is given here - to facilitate the modifications to Your Public Environment
parsePublicLocations(){

# THE ORDER THE ABOVE TAKE EFFECT :

# 1) REMOVE from public path
# 2) ADD to public path
# 3) REPLACE in public path

echo "Parse server path from disk path .."
echo "Parse server path from disk path .." >&2

echo
echo "diskFolder ${diskFolder}"
# The diskFolder is the second column of PIPE_hubbing.sh or PIPE_hubbingSymbolic.txt

serverFolderTEMP=${diskFolder}

# 1) REMOVE from public path

if [ "${REMOVEfromPUBLICFILEPATH}" != "" ];then
    
# Replace all / with \/  to be readable for SED parser
removeParseForSed=$( echo ${REMOVEfromPUBLICFILEPATH} | sed 's/\//\\\//g' )
# Remove from the beginning (^) of the disk path
removed_from_path=$( echo ${serverFolderTEMP} | sed 's/^'${removeParseForSed}'//' )
# Update the serverFolderTEMP
serverFolderTEMP=${removed_from_path}

fi

# 2) ADD to public path

if [ "${ADDtoPUBLICFILEPATH}" != "" ];then

# Replace all / with \/  to be readable for SED parser
addParseForSed=$( echo ${ADDtoPUBLICFILEPATH} | sed 's/\//\\\//g' )
# Add to the beginning (^) of the disk path
added_to_path=$( echo ${serverFolderTEMP} | sed 's/^/'${addParseForSed}'/' )
# Update the serverFolderTEMP
serverFolderTEMP=${added_to_path}

fi

# 3) REPLACE in public path

if [ "${tobeREPLACEDinPUBLICFILEPATH}" != "" ];then
    
# Replace all / with \/  to be readable for SED parser
toBeParseForSed=$( echo ${tobeREPLACEDinPUBLICFILEPATH} | sed 's/\//\\\//g' )
withThisParseForSed=$( echo ${REPLACEwithThisInPUBLICFILEPATH} | sed 's/\//\\\//g' )
# Replace within the disk path
replaced_path=$( echo ${serverFolderTEMP} | sed 's/'${toBeParseForSed}'/'${withThisParseForSed}'/' )
# Update the serverFolderTEMP
serverFolderTEMP=${replaced_path}

# Note !! - the (3) is vulnerable to repeating of the same subfolder structure.
# It only will find and replace the LEFTmost of the occurrences of the file path it searches for.

fi

# This is the value which is returned to the script :

serverFolder="${serverFolderTEMP}"

echo
echo "serverFolder ${serverFolder}"
    
}
