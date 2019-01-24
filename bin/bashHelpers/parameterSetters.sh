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


setMparameter(){
   
mParameter=""

if [ "${CAPITAL_M}" -eq 0 ] ; then
    mParameter="-m ${LOWERCASE_M}"
else
    mParameter="-M ${CAPITAL_M}"
fi 
    
}


setParameters(){

#----------------------------------------------
# Listing current limitations, exiting if needed :

if [ "${LOWERCASE_M}" -ne 0 ] && [ "${CAPITAL_M}" -ne 0 ];
then 
    printThis="Bowth -m and -M parameters cannot be set at the same time\nEXITING"
    printToLogFile
   exit 1
fi

#----------------------------------------------

if [ "${LOWERCASE_V}" -ne -1 ] && [ "${bowtie1MismatchBehavior}" != "" ]
then
    printThis="Bowtie1 does not allow setting -v with any other mismatch-reporting altering parameters ( --seedmms --seedlen --maqerr ) \nUse only -v, or (any) combination of --seedmms --seedlen --maqerr\nEXITING"
    printToLogFile
  exit 1
fi

if [ "${bowtie1MismatchBehavior}" != "" ]
then
    otherBowtie1Parameters="${otherBowtie1Parameters} ${bowtie1MismatchBehavior}"
fi

if [ "${LOWERCASE_V}" -ne -1 ]
then
    otherBowtie1Parameters="${otherBowtie1Parameters} -v ${LOWERCASE_V}"
fi


otherBowtie2Parameters="${otherBowtie2Parameters} ${bowtie2MismatchBehavior}"

#----------------------------------------------
#Setting the m and M parameters..

if [ "${LOWERCASE_M}" -ne 0 ] ;
then
   CAPITAL_M=0 
fi

if [ "${CAPITAL_M}" -ne 0 ];
then
   LOWERCASE_M=0
fi

if [ "${LOWERCASE_M}" -eq 0 ] && [ "${CAPITAL_M}" -eq 0 ];
then
    LOWERCASE_M=2
fi
#----------------------------------------------
#Setting the dependencies of FOOTPRINT and PEAKCALL parameters in rerun..

if [ "${ONLY_PEAK}" -ne 0 ] ;
then
   PEAKCALL=1
fi

if [ "${ONLY_FP_AND_PEAK}" -ne 0 ];
then
   PEAKCALL=1
   FOOTPRINT=1
fi


#----------------------------------------------

# Setting WINDOWING PARAMETERS to be compatible - if they weren't when called

if [ "${DOWINDOWING}" -eq 1 ] ;
then
   
    # If we have only one of those (or neither) set :
    
    # This first if means "default increment" - i.e. 1/10 increment
    if [ "${WINDOWTRACK}" -ne 0 ]
    then
    # Window width has to be dividable by 10.
    WINDOWTRACK=$(($((${WINDOWTRACK}/10))*10))
    # Increment has to divide evenly to track.
    WINDOWINCR=$((${WINDOWTRACK}/10))
    
    # This second if means "default window width" - i.e. 300 bases
    elif [ "${WINDOWTRACK}" -ne 0 ]
    then
    WINDOWTRACK=300
    WINDOWINCR=$(($((${WINDOWINCR}/10))*10))
    
    else
    # Default window
    WINDOWTRACK=300
    WINDOWINCR=30
    
    fi
    
fi

#------------------------------------------------
# Custom adapter sequences..

if [ "${ADA31}" != "no"  ] || [ "${ADA32}" != "no" ] || [ "${ADA51}" != "no" ] || [ "${ADA52}" != "no" ] 
    then
    CUSTOMAD=1
fi

}

setBowtieReadList(){
    
if [ "${SINGLE_END}" -eq 0 ] ; then
    bowtieReadList="-1 ./READ1.fastq -2 ./READ2.fastq"
else
    if [ "${BOWTIE}" -eq 2 ] ; then
        bowtieReadList="-U ./READ1.fastq"
    else
        bowtieReadList="./READ1.fastq"
    fi  
        
fi

}

