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

fetchFastq(){
    
    printThis="Fetching fastq files  ..."
    printNewChapterToLogFile
    
    doQuotaTesting
    
    printThis="${fileList1[$i]}"
    printToLogFile
    testedFile="${fileList1[$i]}"
    doInputFileTesting
    
    echo "original size (R1)"
    ls -lht ${fileList1[$i]} | cut -d " " -f 1,2,3,4 --complement
    
    if [ "${GZIP}" -eq 0 ] ; then
    
    echo "cp ${fileList1[$i]} ./READ1.fastq"
    cp "${fileList1[$i]}" ./READ1.fastq
    
    else
    
    echo "cp ${fileList1[$i]} ./READ1.fastq.gz"
    cp "${fileList1[$i]}" ./READ1.fastq.gz
    echo "original size (R1 copied)"
    ls -lht READ1.fastq.gz | cut -d " " -f 1,2,3,4 --complement
    
    gzip -d READ1.fastq.gz
    
    fi
    
    testedFile="READ1.fastq"
    doTempFileTesting
    
    #echo "ln -fs ${fileList1[$i]} READ1.fastq"
    #ln -fs "${fileList1[$i]}" READ1.fastq
    
    if [ "${SINGLE_END}" -eq 0 ] ; then
    
    printThis="${fileList2[$i]}"
    printToLogFile
    testedFile="${fileList2[$i]}"
    doInputFileTesting
    
    echo "original size (R2)"
    ls -lht ${fileList2[$i]} | cut -d " " -f 1,2,3,4 --complement

    if [ "${GZIP}" -eq 0 ] ; then
    
    echo "cp ${fileList2[$i]} ./READ2.fastq"
    cp "${fileList2[$i]}" ./READ2.fastq
    
    else
    
    echo "cp ${fileList2[$i]} ./READ2.fastq.gz"
    cp "${fileList2[$i]}" ./READ2.fastq.gz
    echo "original size (R2 copied)"
    ls -lht READ2.fastq.gz | cut -d " " -f 1,2,3,4 --complement
    
    gzip -d READ2.fastq.gz
    
    fi
    
    testedFile="READ2.fastq"
    doTempFileTesting
    #echo "ln -fs ${fileList2[$i]} READ2.fastq"
    #ln -fs "${fileList2[$i]}" READ2.fastq
    
    fi

    echo "Fetched fastqs :"
    ls -lh | grep fastq | cut -d " " -f 1,2,3,4 --complement
    
    doQuotaTesting
    
}

fetchFastqMultilane(){
    
    #folderList=($( cut -f 1 ./PIPE_fastqPaths.txt ))
    #fileList1=($( cut -f 2 ./PIPE_fastqPaths.txt ))
    #fileList2=($( cut -f 3 ./PIPE_fastqPaths.txt ))
    
    printThis="Fetching fastq files (we have $LANES lanes)..."
    printNewChapterToLogFile

    doQuotaTesting
    
    echo "Read1 - generating combined fastq.."
    printThis="${fileList1[$i]}"
    printToLogFile
    

    # One lane at a time.. catenating files !
    rm -f ./READ1.fastq 
    allLanes=${fileList1[$i]}

    for l in $( seq 1 ${LANES[@]} ); do
        echo ""
        echo "Lane no $l .."
        currentLane=$( echo ${allLanes} | sed s'/,.*$//' )
        echo "Current lane : ${currentLane}"
        testedFile="${currentLane}"
        doInputFileTesting
        
        echo "original size (current lane R1)"
        ls -lht ${currentLane} | cut -d " " -f 1,2,3,4 --complement
        if [ "${GZIP}" -eq 0 ] ; then
           cat "${currentLane}" >> ./READ1.fastq
        else
           cp "${currentLane}" ./TEMP.fastq.gz
           echo "original size (current lane R1 copied)"
           ls -lht TEMP.fastq.gz | cut -d " " -f 1,2,3,4 --complement
           gzip -d TEMP.fastq.gz
           echo "unpacked size (current lane R1 unpacked)"
           ls -lht TEMP.fastq | cut -d " " -f 1,2,3,4 --complement
           cat TEMP.fastq >> ./READ1.fastq
           rm -f TEMP.fastq*
        fi
        echo "current size of combined file R1"
        ls -lht READ1.fastq | cut -d " " -f 1,2,3,4 --complement
        
        # Saving rest and looping to next round..
        removeThis=$( echo ${currentLane} | sed 's/\//\\\//g' )
        restOfLanes=$( echo ${allLanes} | sed s'/'${removeThis}',//' )
        echo "Rest of lanes (still to be added to the file) : ${restOfLanes}"
        allLanes=${restOfLanes}  
    done
    # Removing empty lines
    grep -Pv "^$" READ1.fastq >  temp.fastq
    mv -f temp.fastq READ1.fastq
    
    testedFile="READ1.fastq"
    doTempFileTesting
    
    doQuotaTesting
    
    if [ "${SINGLE_END}" -eq 0 ] ; then
    
    echo ""
    echo "Read2 - generating combined fastq.."
    printThis="${fileList2[$i]}"
    printToLogFile

    rm -f ./READ2.fastq 
    allLanes=${fileList2[$i]}
    for l in $( seq 1 ${LANES[@]} ); do
        echo ""
        echo "Lane no $l .."
        currentLane=$( echo ${allLanes} | sed s'/,.*$//' )
        echo "Current lane : ${currentLane}"
        testedFile="${currentLane}"
        doInputFileTesting
        
        echo "original size (current lane R2)"
        ls -lht ${currentLane} | cut -d " " -f 1,2,3,4 --complement
        if [ "${GZIP}" -eq 0 ] ; then
           cat "${currentLane}" >> ./READ2.fastq
        else
           cp "${currentLane}" ./TEMP.fastq.gz
           echo "original size (current lane R2 copied)"
           ls -lht TEMP.fastq.gz | cut -d " " -f 1,2,3,4 --complement
           gzip -d TEMP.fastq.gz
           echo "unpacked size (current lane R2 unpacked)"
           ls -lht TEMP.fastq | cut -d " " -f 1,2,3,4 --complement
           cat TEMP.fastq >> ./READ2.fastq
           rm -f TEMP.fastq*
        fi
        echo "current size of combined file R2"
        ls -lht READ2.fastq | cut -d " " -f 1,2,3,4 --complement
        
        # Saving rest and looping to next round..
        removeThis=$( echo ${currentLane} | sed 's/\//\\\//g' )
        restOfLanes=$( echo ${allLanes} | sed s'/'${removeThis}',//' )
        echo "Rest of lanes (still to be added to the file) : ${restOfLanes}"
        allLanes=${restOfLanes}   
    done
    # Removing empty lines
    grep -Pv "^$" READ2.fastq >  temp.fastq
    mv -f temp.fastq READ2.fastq
    
    testedFile="READ2.fastq"
    doTempFileTesting
    
    doQuotaTesting
    
    fi

    echo "Generated merged fastqs :"
    ls -lh | grep fastq | cut -d " " -f 1,2,3,4 --complement
    
}

generateFastq(){
    #Needs this to be set before call
    #bamPath=/path/to/file.bam
    
    printThis="Fetching files (and converting to fastqs) ..."
    printNewChapterToLogFile
    
    doQuotaTesting
    
    printThis="${bamPath}"
    printToLogFile

    samtools sort -n "${bamPath}" SortedByName
    if [ ! -r "SortedByName.bam" ] || [ ! -e "SortedByName.bam" ] || [ ! -f "SortedByName.bam" ] ; then
      echo "Temporary file not found or empty file : SortedByName.bam" >&2
      echo "EXITING!!" >&2
      exit 1
    fi
 
    samtools view -o ForBam2Fastq.sam SortedByName.bam
    if [ ! -r "ForBam2Fastq.sam" ] || [ ! -e "ForBam2Fastq.sam" ] || [ ! -f "ForBam2Fastq.sam" ] ; then
      echo "Temporary file not found or empty file : ForBam2Fastq.sam" >&2
      echo "EXITING!!" >&2
      exit 1
    fi
    
    rm -f SortedByName.bam
    printThis=" perl ${PerlHelpersPath}/sam2fastq.pl < ForBam2Fastq.sam > sam2fastq.log"
    printToLogFile
    perl ${PerlHelpersPath}/sam2fastq.pl < ForBam2Fastq.sam > sam2fastq.log
    
    rm -f ForBam2Fastq.sam
    
    doQuotaTesting

    echo "Generated fastqs :"
    ls -lh | cut -d " " -f 1,2,3,4 --complement
    echo "In folder :"
    pwd
}

