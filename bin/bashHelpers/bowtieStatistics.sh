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

doStatisticsFile(){
    
    # This subroutine is copied from afterBowtieMapping.sh on 110714 14:28
    
    #Demands coordinate-wise sorted and indexed bam file
    
    # NEEDS THIS TO BE SET BEFORE CALL :
    # thisSortedBAMfile="" - file name WITHOUT ./ as it will enter SED 
    # Assumes also that c.o. file is also INDEXED already. 
    
    echo "STATISTICS file generation for file : ${thisSortedBAMfile}" >&2
    
    tempName=$(echo ${thisSortedBAMfile} | sed 's/_Sorted.bam//')
    logName="${tempName}_statistics.log"
    
    echo "Total amount of reads from samtools view -c : " > ${logName}       
    samtools view -c ${thisSortedBAMfile} >> ${logName}              
    echo "" >> ${logName}
    
    echo FLAGSTAT >> ${logName}
    samtools flagstat ${thisSortedBAMfile} >> ${logName}
    echo "" >> ${logName}
    echo IDXSTATS >> ${logName}
    samtools idxstats ${thisSortedBAMfile} >> ${logName}
}


