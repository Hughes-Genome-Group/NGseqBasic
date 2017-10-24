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

listFolder(){
    
        echo "-----------------------------------------------"
        echo "Now finished in folder"
        pwd
        echo
        echo "Generated files and folders :"
        ls -lh | cut -d " " -f 1,2,3,4 --complement
        echo
        echo "Bigwig files in BigWigs folder :"
        ls -lh BigWigs | cut -d " " -f 1,2,3,4 --complement
        echo
        
        if [ -d BedGraphs_and_Wigs ]; then
        
        echo "Bedgraph and wig files in BedGraphs_and_Wigs folder :"
        ls -lh BedGraphs_and_Wigs | cut -d " " -f 1,2,3,4 --complement
        echo
        
        fi
        
        echo "Samtools statistics files in SamTools_statistics folder :"
        ls -lh SamTools_statistics | cut -d " " -f 1,2,3,4 --complement
        echo "-----------------------------------------------"
        echo

}
    
cleanUpFolder() {
    
        # Cleaning up in the folder !
        mkdir SamTools_statistics BedGraphs_and_Wigs BigWigs
        mv -f *statistics.log SamTools_statistics/ 2> /dev/null
        mv -f *.bw BigWigs/ 2> /dev/null
        mv -f *.bdg BedGraphs_and_Wigs/ 2> /dev/null
        mv -f *.wig BedGraphs_and_Wigs/ 2> /dev/null
        
        rm -f extendedFrags.fastq
        rm -f flash.hist
        rm -f flash.histogram
        
        if [ "${SAVE_BDG_AND_WIG}" -eq 0 ]; then
            rm -rf BedGraphs_and_Wigs
        fi
    
}

deleteFastqLocal(){
   
    #If problems in bowtie - producing empty or non-existent output file, aborting script :
        testedFile="bowtie_out.bam"
        doTempFileTesting
  
    printThis="Deleting the fastq SUBFOLDER files.."
    printNewChapterToLogFile
    
    rm -f READ1.fastq
    rm -f READ2.fastq

echo ""

}

deleteFastq(){
   
    #If problems in bowtie - producing empty or non-existent output file, aborting script :
    
    
    for k in $( seq 0 $((${#GENOMEARRAY[@]} - 1)) ); do
    #for k in "${GENOMEARRAY[@]}"; do
        genomeName=${GENOMEARRAY[$k]}
        testedFile="${genomeName}/bowtie_out.bam"
        doTempFileTesting
    done
  
    printThis="Deleting the fastq MAIN FOLDER files.."
    printNewChapterToLogFile
    
    rm -f READ1.fastq
    rm -f READ2.fastq

echo ""

}

