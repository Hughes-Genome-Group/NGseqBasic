#!/bin/bash

##########################################################################
# Copyright 2016, Jelena Telenius (jelena.telenius@imm.ox.ac.uk)         #
#                                                                        #
# This file is part of NGseqBasic .                                      #
#                                                                        #
# NGseqBasic is free software: you can redistribute it and/or modify     #
# it under the terms of the MIT license.
#
#
#                                                                        #
# NGseqBasic is distributed in the hope that it will be useful,          #
# but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
# MIT license for more details.
#                                                                        #
# You should have received a copy of the MIT license
# along with NGseqBasic.  
##########################################################################

doInputFileTesting(){
    
    # NEEDS THIS TO BE SET BEFORE CALL :
    # testedFile=""    
    
    if [ ! -r "${testedFile}" ] || [ ! -e "${testedFile}" ] || [ ! -f "${testedFile}" ] || [ ! -s "${testedFile}" ]; then
      echo "Input file not found or empty file : ${testedFile}" >> "/dev/stderr"
      echo "EXITING!!" >> "/dev/stderr"
      exit 1
    fi
}

doTempFileTesting(){
    
    # NEEDS THIS TO BE SET BEFORE CALL :
    # testedFile=""    
    
    if [ ! -r "${testedFile}" ] || [ ! -e "${testedFile}" ] || [ ! -f "${testedFile}" ] || [ ! -s "${testedFile}" ]; then
      echo "Temporary file not found or empty file : ${testedFile}" >&2
      echo "EXITING!!" >&2
      exit 1
    fi
}

doTempFileInfo(){
    
    # NEEDS THIS TO BE SET BEFORE CALL :
    # testedFile=""    
    
    if [ ! -r "${testedFile}" ] || [ ! -e "${testedFile}" ] || [ ! -f "${testedFile}" ] || [ ! -s "${testedFile}" ]; then
      echo "Temporary file not found or empty file : ${testedFile}" >&2
      echo "Continuing, but may create some NONSENSE data.." >&2
      #echo "EXITING!!" >&2
      #exit 1
    fi
}

doTempFileFYI(){
    
    # NEEDS THIS TO BE SET BEFORE CALL :
    # testedFile=""    
    
    if [ ! -r "${testedFile}" ] || [ ! -e "${testedFile}" ] || [ ! -f "${testedFile}" ] || [ ! -s "${testedFile}" ]; then
      echo "Probably harmless : file not found or empty file : ${testedFile}" >&2
      
      #echo "EXITING!!" >&2
      #exit 1
    fi
}

doQuotaTesting(){
    
    echo
    if [ $( pwd | grep -c '^/t1-data1/' ) -ne 0 ]; then
        getquota-t1data1 2>>/dev/null
    else
        getquota 2>>/dev/null
    fi
        
    echo
    echo "Disk usage for THIS RUN - at the moment :"
    du -sh ${dirForQuotaAsking}
    echo
    
}

