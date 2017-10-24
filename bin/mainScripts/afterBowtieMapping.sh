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


# This code was updated to current DNasePipe version 241216

doBigWigGeneration(){
    echo BIGWIG GENERATION >&2

    # Note that the "whole read pair" bigwigs and bedgraphs contain SINGLE INTERVAL for each properly paired READ PAIR.
    # The footprint bigwigs and bedgraphs contain 1bp from the LEFT and RIGHT ends of the read (respectively)
    
    # Red and green graphs (all reads and filtered reads)

    thisBDGfile="filtered_pileup.bdg"
    DoRounding=0
    doBigWigging
    
    thisBDGfile="bowtie_out_mapped_pileup.bdg"
    DoRounding=0
    doBigWigging
    
    if [ "${WINTRACK}" -ne 0 ] ; then
        
    thisBDGfile="filtered_window.bdg"
    DoRounding=0
    doBigWigging
    
    thisBDGfile="bowtie_out_mapped_window.bdg"
    DoRounding=0
    doBigWigging
    
fi

    thisIsWhereIam=$( pwd )
    printThis="Starting to sort filtered_pileup.bdg BIG time - will save temporary files in ${thisIsWhereIam}"
    printToLogFile
    
    if [ "${WINTRACK}" -eq 0 ] ; then
    
    cut -f 4 filtered_pileup.bdg |  awk '{s=(($1)/10)} {print s}' > tempCol4.txt
    cut -f 1,2,3 filtered_pileup.bdg | paste - tempCol4.txt > intoSorting.txt
    
    else
    
    cut -f 4 filtered_window.bdg |  awk '{s=(($1)/100)} {print s}' > tempCol4.txt
    cut -f 1,2,3 filtered_window.bdg | paste - tempCol4.txt > intoSorting.txt

    fi
   
    sortParams='-k1,1 -k2,2n'
    sortIn1E6bunches
  # needs these to be set :
  # thisIsWhereIam=$( pwd )
  # sortParams="-k1,1 -k2,2n"  or sortParams="-n" etc
  # input in intoSorting.txt
  # outputs TEMPsortedMerged.txt
    
    # If all went well, we delete original file. If not, we complain here (but will not die).
   sortResultInfo
   rm -f intoSorting.txt
    
   mv -f TEMPsortedMerged.txt filtered_pileup_scaled.bdg   
    
    
    thisBDGfile="filtered_pileup_scaled.bdg"
    DoRounding=1
    doBigWigging
    rm -f filtered_pileup_scaled.bdg tempCol4.txt
    
    
    
    # Unpaired mapped
    
    thisBDGfile="singleEnd_READ1_filtered_pileup.bdg"
    DoRounding=0
    doBigWigging
    
    thisBDGfile="singleEnd_READ2_filtered_pileup.bdg"
    DoRounding=0
    doBigWigging

    thisBDGfile="singleEnd_READ1_pileup.bdg"
    DoRounding=0
    doBigWigging
    
    thisBDGfile="singleEnd_READ2_pileup.bdg"
    DoRounding=0
    doBigWigging

    
    if [ "${FOOTPRINT}" -eq 1 ] ; then
        
    thisBDGfile="LEFT_pileup.bdg"
    DoRounding=0
    doBigWigging

    thisBDGfile="RIGHT_pileup.bdg"
    DoRounding=0
    doBigWigging
    
    thisBDGfile="LEFT_window.bdg"
    DoRounding=1
    doBigWigging
    rm -f LEFT_window.bdg

    thisBDGfile="RIGHT_window.bdg"
    DoRounding=1
    doBigWigging
    rm -f RIGHT_window.bdg
    
    fi


    if [ "${PEAKCALL}" -eq 1 ] ; then
    
    thisBDGfile="READ_peakCall.bdg"
    DoRounding=0
    doBigWigging
    
    thisBDGfile="READ_peakCall_quant.bdg"
    DoRounding=0
    doBigWigging
    
    echo -e "chr\tstart\tstop\treadDepth\tpeakWidth" > heading.txt
    cat READ_peakCall_quant.bdg | cut -f 2,3 | awk '{ print $2-$1 }' > TEMPwidth.txt
    paste READ_peakCall_quant.bdg TEMPwidth.txt > TEMP_intoMIG.txt
    cat heading.txt TEMP_intoMIG.txt > forMigging.txt
    rm -f TEMPwidth.txt TEMP_intoMIG.txt
    
    perl ${PerlHelpersPath}/data2gff.pl --data forMigging.txt --output PeakCall.gff
    rm -f heading.txt forMigging.txt
    
    fi
    
    echo
    echo "Generated files :"
    ls -lh | cut -d " " -f 1,2,3,4 --complement
    echo
    echo
    
}

doStatisticsFile(){
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

doBigWigging(){
    
    # NEEDS THIS TO BE SET BEFORE CALL :
    # thisBDGfile=""
    # DoRounding="" (0 or 1)
    
    if [ -r "${thisBDGfile}" ] ; then
    if [ -s "${thisBDGfile}" ] ; then
    bigwigName=$(echo ${thisBDGfile} | sed 's/\..*//')
    
    if [ "${DoRounding}" -eq 0 ]; then
    
    bedGraphToBigWig ${thisBDGfile} ${genomeBuild} "${bigwigName}.bw"

    
    #Temporary file to be bigwigged - the awk command here drops the resolution of the bigwig to be 2 digits..
    else
    
    #cut -f 4 ${thisBDGfile} | grep -Pv "^$" | awk '{printf("%.2f",$1);  printf("\n");}' > tempCol4.txt
    cut -f 4 ${thisBDGfile} | awk '{printf("%.2f",$1);  printf("\n");}' > tempCol4.txt
    
    #tail tempCol4.txt
    
    cut -f 1,2,3 ${thisBDGfile} | paste - tempCol4.txt > toBigWig.bdg
    
    #head toBigWig.bdg
    #tail toBigWig.bdg
    
    bedGraphToBigWig toBigWig.bdg ${genomeBuild} "${bigwigName}.bw"
    rm -f tempCol4.txt toBigWig.bdg
    
    fi
    
    #ls -lhs| grep ${bigwigName} >&2
    
    else
    echo "File ${thisBDGfile} empty file, no bigwig generated !"  >&2
    fi
    else
    echo "File ${thisBDGfile} missing, no bigwig generated !"    >&2 
    fi
}


doBigBedding(){
    
    # NEEDS THIS TO BE SET BEFORE CALL :
    # thisBEDfile=""
    
    if [ -r "${thisBEDfile}" ] ; then
    if [ -s "${thisBEDfile}" ] ; then
    bigbedName=$(echo ${thisBEDfile} | sed 's/\..*//')
    
    bedToBigBed ${thisBEDfile} ${genomeBuild} "${bigbedName}.bb"
    
    else
    echo "File ${thisBEDfile} empty file, no bigbed generated !"  >&2
    fi
    else
    echo "File ${thisBEDfile} missing, no bigbed generated !"    >&2 
    fi
}

sortIn1E6bunches(){
  # needs these to be set :
  # thisIsWhereIam=$( pwd )
  # sortParams="-k1,1 -k2,2n"  or sortParams="-n" etc
  # input in intoSorting.txt
  # outputs TEMPsortedMerged.txt
  
  rm -f preSorting.txt
  cp intoSorting.txt preSorting.txt
  
   # Here we make the files in 1 000 000 reads bunches
    
   head -n  1000000 preSorting.txt > forSortingfile1.txt
   tail -n +1000001 preSorting.txt > temp.txt
   mv -f temp.txt preSorting.txt
   
   TEMPPcounter=1
   while [ -s preSorting.txt ]
   do
        TEMPPcounter=$(( ${TEMPPcounter}+1 ))
        head -n  1000000 preSorting.txt > forSortingfile${TEMPPcounter}.txt
        tail -n +1000001 preSorting.txt > temp.txt
        mv -f temp.txt preSorting.txt
        
   done
   rm -f preSorting.txt
   
   echo "made this many files for sorting (each upto 1 million reads) :"  >&2
   # ls -lht | grep forSortingfile >&2
   ls forSortingfile* >&2
   
   for file in forSortingfile*
   do
       newwwName=$( echo ${file} | sed 's/forSortingfile/preSortedfile/' )
       # echo "sort -S ${memoryMegas}M ${sortParams} -T ${thisIsWhereIam} ${file} > ${newwwName}" >&2
       sort -S ${memoryMegas}M ${sortParams} -T ${thisIsWhereIam} ${file} > ${newwwName}
       rm -f ${file}
       
   done
   
   sort -m -S ${memoryMegas}M ${sortParams} -T ${thisIsWhereIam} preSortedfile* > TEMPsortedMerged.txt
   rm -f preSortedfile*
   
}

sortResultTester(){
    
     # Check if all went well.. 

    countBefore=$(($( cat intoSorting.txt | grep -c "" )))
    countAfter=$(($( cat TEMPsortedMerged.txt | grep -c "" )))
    
    if [ "${countBefore}" -ne "${countAfter}" ]
    then
    
    echo "BIG TIME Sorting FAILED. " >&2
    echo "Original file had ${countBefore} data lines - sorted file had only ${countAfter} lines." >&2
    echo "EXITING!!" >&2
    exit 1
    
    else

    echo "BIG TIME Sorting SUCCEEDED ! " >&2

    fi   
    
}

sortResultInfo(){
    
    # Check if all went well.. 

    countBefore=$(($( cat intoSorting.txt | grep -c "" )))
    countAfter=$(($( cat TEMPsortedMerged.txt | grep -c "" )))
    
    if [ "${countBefore}" -ne "${countAfter}" ]
    then
    
    echo "BIG TIME Sorting FAILED." >&2
    echo "Original file had ${countBefore} data lines - sorted file had only ${countAfter} lines." >&2
    echo "CONTINUING - but may produce nonsense data!!" >&2
    
    else

    echo "BIG TIME Sorting SUCCEEDED ! " >&2

    fi
}

doWindowing(){
    # NEEDS THESE TO BE SET BEFORE CALL :
    # WindowInput="" the name of the file)
    # WindowOutputName="" (identifier of the data file - to be outputted in output too)
    # WINDOW=""
    # INCREMENT=""
    # genomeBuild=""    
    
# Do only if there are no SCAFFOLDS in the genome..

if [ ! "${genomeName}" = "xenLae9" ] ; then
    
    testedFile="${WindowInput}"
    doTempFileInfo
    
    testedFile="${genomeBuild}"
    doTempFileInfo
    
#--------------------------------------------------------------------------------------------------------------------------------
printThis="${WindowInput} -running smoothing over the data set ${WINDOW} bp window, ${INCREMENT} bp overlap"
printToLogFile

#--------------------------------------------------------------------------------------------------------------------------------

# This is heavy duty stuff - doing it line by line..

#bedtools makewindows -w $((${WINDOW}*4)) -s $((${WINDOW}*2)) ${WindowType} ${WindowRegions} > windowOver3.bed
#bedtools coverage -counts -a ${WindowInput} -b windowOver3.bed > coverageOver3.txt

cut -f 1-3 ${WindowInput} > TEMP_input.bed

rm -rf TEMPwindow
mkdir TEMPwindow
cd TEMPwindow
    cat -n ${genomeBuild} | awk '{printf "%s\t%s\t%s\n", $2,$3,$4 >$1".seq"; close($1".seq")}' 
cd ..

regioncount=$( ls -1 TEMPwindow | grep -c "" )

rm -f windowOver3.bed coverageOver3.txt
echo -n "" > coverageOver3.txt

echo
echo -n "Making windows.. Region "

regioncounter=1
for windowFile in TEMPwindow/*
do
    echo -n "${regioncounter}/${regioncount}.. "
    regioncounter=$(( ${regioncounter}+1 ))
    
    # The true windowing is done sam2bw-style in perl script :
    
    perl ${PerlHelpersPath}/windowingScript.pl -window ${WINDOW} -inc ${INCREMENT} -genome ${windowFile} -bed TEMP_input.bed -name TEMP_perl_output
    cat TEMP_perl_output.windowed.bed >>  coverageOver3.txt
    rm -f TEMP_perl_output.windowed.bed

done
echo

rm -rf TEMPwindow TEMP_input.bed

testedFile="coverageOver3.txt"
doTempFileInfo

thisIsWhereIam=$( pwd )
printThis="Starting to sort coverageOver3.txt BIG time - will save temporary files in ${thisIsWhereIam}"
printToLogFile

mv -f coverageOver3.txt intoSorting.txt
    sortParams='-k1,1 -k2,2n'
    sortIn1E6bunches
  # needs these to be set :
  # thisIsWhereIam=$( pwd )
  # sortParams="-k1,1 -k2,2n"  or sortParams="-n" etc
  # input in intoSorting.txt
  # outputs TEMPsortedMerged.txt
    
    # If all went well, we delete original file. If not, we complain here (but will not die).
   sortResultInfo
   rm -f intoSorting.txt

   mv -f TEMPsortedMerged.txt beforeClip.bed

bedClip beforeClip.bed ${genomeBuild} beforeCompress.bdg
rm -f beforeClip.bed
bedGraphPack beforeCompress.bdg "${WindowOutputName}"
rm -f beforeCompress.bdg

testedFile="${WindowOutputName}"
doTempFileInfo

#    echo "--------------------------------------------------"
#    ls -lh | cut -d " " -f 1,2,3,4 --complement  >&2


# End if for "if scaffold genome"
else
    echo
    echo "Skipping windowing for track ${WindowInput} - genome ${genomeName} contains SCAFFOLDS, and those are not supported by the windowing script."
    echo
fi

}

unpairedBdg(){
    # NEEDS THIS TO BE SET BEFORE CALL :
    # unpairedSorted=""
    # Assumes that the file is already indexed, as well.
    
    thisSortedBAMfile=${unpairedSorted}
    doStatisticsFile
    
    baseName=$(echo ${unpairedSorted} | sed 's/Sorted\.bam//')
    
    #UNFILTERED unpaired
    echo UNFILTERED unpaired GENERATION >&2
    outputName="${baseName}pileup.bdg"
    thisIsWhereIam=$( pwd )
    printThis="Starting to sort ${unpairedSorted} BIG time - will save temporary files in ${thisIsWhereIam}"
    printToLogFile
    bedtools bamtobed -i ${unpairedSorted} | cut -f 1,2,3 > intoSorting.txt
    
    sortParams='-k1,1 -k2,2n'
    sortIn1E6bunches
  # needs these to be set :
  # thisIsWhereIam=$( pwd )
  # sortParams="-k1,1 -k2,2n"  or sortParams="-n" etc
  # input in intoSorting.txt
  # outputs TEMPsortedMerged.txt
    
    # If all went well, we delete original file. If not, we complain here (but will not die).
   sortResultInfo
   rm -f intoSorting.txt
   
   cat TEMPsortedMerged.txt | bedtools genomecov -bg -i stdin -g ${genomeBuild} > ${outputName}
   rm -f TEMPsortedMerged.txt
    
    #fi
    
    #FILTERED unpaired
    echo FILTERED unpaired GENERATION >&2
    samtools rmdup -s ${unpairedSorted} nodup_unpaired.bam
    
    if [ "$(($( samtools view -c nodup_unpaired.bam )))" -eq 0 ]; then
      echo "No unpaired reads in filtered unpaired file !" >&2
      echo "No filtered_unpaired.bam file or its visualisation generated"
      rm -f nodup_unpaired.bam
    else
    
    #echo "--------------------------------------------------"
    #ls -lh | cut -d " " -f 1,2,3,4 --complement  >&2
    
    
    # HERE TO BE FLAGGED - IF THESE ARE WISHED TO BE SAVED ! (saved when flag value 11 or 1)
    if [ "${saveUnpaired}" -eq 0 ] || [ "${saveUnpaired}" -eq 10 ]; then
        rm -f ${unpairedSorted}*
    else
        samtools index ${unpairedSorted}
    fi
    
    samtools sort -n nodup_unpaired.bam nodup_unmapped_Sorted
    
    if [ ! -r "./nodup_unmapped_Sorted.bam" ] || [ ! -e "./nodup_unmapped_Sorted.bam" ] || [ ! -f "./nodup_unmapped_Sorted.bam" ]; then
      echo "Temporary file not found or empty file : nodup_unmapped_Sorted.bam" >&2
      echo "No filtered_unpaired.bam file or its visualisation generated"
    else
        
    #echo "--------------------------------------------------"
    #ls -lh | cut -d " " -f 1,2,3,4 --complement  >&2
    rm -f nodup_unpaired.bam
    
    # HERE TO BE REMOVED (IF_CLAUSE) WHEN SUPPORTING OTHER PLOIDY TRACKS TOO !
    if [ "${genomeName}" = "mm9" ] || [ "${genomeName}" = "mm10" ] || [ "${genomeName}" = "hg18" ] || [ "${genomeName}" = "hg19" ] ; then 
    bedtools intersect -v -sorted -abam nodup_unmapped_Sorted.bam -b "${ploidyPath}" > filtered_unpaired.bam
    else
    mv nodup_unmapped_Sorted.bam filtered_unpaired.bam    
    fi
    
    #echo "--------------------------------------------------"
    #ls -lh | cut -d " " -f 1,2,3,4 --complement  >&2
    rm -f nodup_unmapped_Sorted.bam
    
    # FILTERED unpaired STATISTICS
    
    outputName2="${baseName}filtered_Sorted"
    if [ ! -r filtered_unpaired.bam ] || [ ! -e filtered_unpaired.bam ] || [ ! -f filtered_unpaired.bam ]; then
      echo "Temporary file not found or empty file : filtered_unpaired.bam" >&2 
      echo "No filtered_unpaired.bam file statistics generated"
    else
        samtools sort filtered_unpaired.bam ${outputName2}
    fi
    
    if [ ! -r "${outputName2}.bam" ] || [ ! -e "${outputName2}.bam" ] || [ ! -f "${outputName2}.bam" ]; then
      echo "Temporary file not found or empty file : ${outputName2}.bam" >&2 
      echo "No filtered_unpaired.bam file statistics generated"
    else
    
    samtools index ${outputName2}.bam
    
    thisSortedBAMfile="${outputName2}.bam"
    doStatisticsFile
    
    # FILTERED unpaired VISUALISATION
    echo FILTERED unpaired VISUALISATION >&2
    
    thisIsWhereIam=$( pwd )
    printThis="Starting to sort filtered_unpaired.bam BIG time - will save temporary files in ${thisIsWhereIam}"
    printToLogFile
    bedtools bamtobed -i filtered_unpaired.bam | cut -f 1,2,3 > intoSorting.txt
    
    sortParams='-k1,1 -k2,2n'
    sortIn1E6bunches
  # needs these to be set :
  # thisIsWhereIam=$( pwd )
  # sortParams="-k1,1 -k2,2n"  or sortParams="-n" etc
  # input in intoSorting.txt
  # outputs TEMPsortedMerged.txt
    
    # If all went well, we delete original file. If not, we complain here (but will not die).
   sortResultInfo
   rm -f intoSorting.txt
   
   mv -f TEMPsortedMerged.txt forCoverages.bed
 
 
     
        outputName3="${baseName}filtered_pileup.bdg"
        bedtools genomecov -bg -i forCoverages.bed -g ${genomeBuild} > ${outputName3}
        
        testedFile="${outputName3}"
        doTempFileInfo
        
    #if [ "${WINTRACK}" -ne 0 ] ; then  
    #    
    #    WindowType="-g"
    #    WindowRegions="${genomeBuild}"
    #    WindowInput="forCoverages.bed"
    #    WindowOutputName="${baseName}filtered_window.bdg"
    #    WINDOW=${WINTRACK}
    #    INCREMENT=${INCTRACK}
    #    doWindowing
    #    
    #fi
    
    rm -f forCoverages.bed
   
    # Removing the name-wise sorted
    #ls -lh | cut -d " " -f 1,2,3,4 --complement  >&2
    rm -f filtered_unpaired.bam*
    
    # HERE TO BE FLAGGED - IF THESE ARE WISHED TO BE SAVED !
    if [ "${saveUnpaired}" -eq 0 ]; then
        rm -f "${outputName2}".bam*
    else
        samtools view -b -F 4 -o "singleEnd_bowtie_filtered_"${baseName}".bam" ${outputName2}".bam"
        samtools index "singleEnd_bowtie_filtered_"${baseName}".bam"
    fi

    #closing all the 3 if of "if unfiltered_unmapped_Sorted.bam" exists
    fi
    fi
    fi
}
doPeakCall(){

# Needs this to be set before calling :
# ThisEnd="LEFT" or "RIGHT"
# MERGEFIRST=1 --> order will be merge, depth, contig
# MERGEFIRST=0 --> order will be depth, merge, contig
# Also needs this file in the running folder : sortedTempFile.bed
# For filtering demands ${MERGE} , ${CONTIG} and ${DEPTH} to be set
# For windowing demands ${WINDOW} to be set
# Generates files ${ThisEnd}_peakCall.bdg, ${ThisEnd}_window.bdg

printThis="${ThisEnd} reads - 1b resolution bedgraph"
printToLogFile
    
    testedFile="sortedTempFile.bed"
    doTempFileInfo
    
    bedtools genomecov -bg -i sortedTempFile.bed -g ${genomeBuild} > "${ThisEnd}_pileup.bdg"
    
if [ "${MERGEFIRST}" -eq 1 ]
then
#--HERE-FOOTPRINT-STYLE-FILTERING-----------------------------------------------------------------------------
    
#--------------------------------------------------------------------------------------------------------------
printThis="Merging.. "
printToLogFile
#--------------------------------------------------------------------------------------------------------------

# First MERGE to really merge regions
# Generating "mergedTempFile.bdg"
cut -f 4 "${ThisEnd}_pileup.bdg" > col4.txt
paste "${ThisEnd}_pileup.bdg" col4.txt > tempFiveCol.bed
bedtools merge -i tempFiveCol.bed -d "${MERGE}" -scores sum > mergedTempFile.bdg
rm -f col4.txt 

#--------------------------------------------------------------------------------------------------------------
printThis="Depth filtering.. "
printToLogFile
#--------------------------------------------------------------------------------------------------------------

# This call needs file depthInFile.bdg
# And produces file peakCall_temp.bdg
# Takes DEPTH as parameter

mv -f mergedTempFile.bdg depthInFile.bdg

awk ' $4 >= '"${DEPTH}"'' depthInFile.bdg > peakCall_temp.bdg

else
#--HERE-TOTAL-READ-DEPTH-STYLE-FILTERING-----------------------------------------------------------------------------
    
#--------------------------------------------------------------------------------------------------------------
printThis="Depth filtering.. "
printToLogFile
#--------------------------------------------------------------------------------------------------------------

awk ' $4 >= '"${DEPTH}"'' "${ThisEnd}_pileup.bdg" > depthOutFile.bdg
    
#--------------------------------------------------------------------------------------------------------------
printThis="Merging.. "
printToLogFile
#--------------------------------------------------------------------------------------------------------------
cut -f 4 depthOutFile.bdg > col4.txt
paste depthOutFile.bdg col4.txt > tempFiveCol.bed
bedtools merge -i tempFiveCol.bed -d "${MERGE}" -scores sum > peakCall_temp_temp.bdg
rm -f col4.txt depthOutFile.bdg
mv -f peakCall_temp_temp.bdg peakCall_temp.bdg
rm -f tempFiveCol.bed peakCall_temp_temp.bdg

fi
    
#--------------------------------------------------------------------------------------------------------------
printThis="Contig filtering.. "
printToLogFile
#--------------------------------------------------------------------------------------------------------------

#Finally, merge complement to get rid of sparse regions(CONTIG) :
complementMerge=$((${CONTIG}-1))

if [ -s peakCall_temp.bdg ]; then

bedtools complement -i peakCall_temp.bdg -g ${genomeBuild} | bedtools merge -d ${complementMerge} -i stdin | bedtools complement -i stdin -g ${genomeBuild}  > mergeContigTempFile.bed
rm -f mergedTempFile.bdg peakCall_temp.bdg

if [ -s "mergeContigTempFile.bed" ]
then
    
sed -i 's/$/\t1/' mergeContigTempFile.bed

mv mergeContigTempFile.bed "${ThisEnd}_peakCall.bdg"

testedFile="${ThisEnd}_peakCall.bdg"
doTempFileInfo

else
    echo "Temporary file mergeContigTempFile.bed empty or not existing (no reads left in file), no pileups generated. Check your --depth --contig and --merge parameters !" >&2

fi

else
    echo "Temporary file peakCall_temp.bdg empty or not existing (no reads left in file), no pileups generated. Check your --depth --contig and --merge parameters !" >&2

fi


#    echo "--------------------------------------------------"
#    ls -lh | cut -d " " -f 1,2,3,4 --complement  >&2
    
}

##################################################################################
#                                                                                #
# ABOVE = the FUNCTION DEFINITIONS (in bash they have to be first in the code)   #
# BELOW = the USAGE information (how to use the code)                            #
# IN THE END = the MAIN CODE itself.                                             #
#                                                                                #
##################################################################################

# Usage - how to use the code :

# 1_bowtie_SubScript_AllGenomes -g mm9
# Run the script one genome build at a time, for one data set at a time. -g parameter is obligatory.

# DOES NOT LOAD MODULES IT NEEDS :
#module load samtools/0.1.19   
#module load bedtools/2.17.0   
#module load ucsctools/1.0

# SHELL SCRIPTS inherit from their parents - (uncommenting) the following line proves that, if you are in doubt !
# module list

# FILES IT REQUIRES :
# bowtie_out.sam - in the folder the script is launched.

# FILES IT DELETES :
# bowtie_out.sam in the run directory.

# FILES IT PRODUCES (in the run directory) :

#Coordinate-wise sorted and indexed bams :
#"""""""""""""""""""""""""""""""""""""""""

#bowtie_out_Sorted.bam
#bowtie_out_Sorted.bam.bai

#filtered_Sorted.bam  (only proper pairs, no duplicates, no ploidy regions)
#filtered_Sorted.bam.bai

#filtered_unmapped_Sorted.bam (mapped, paired in sequencing, not mapped as proper pair - no duplicates, no ploidy regions)
#filtered_unmapped_Sorted.bam.bai

#Log and statistics files for the bams :
#"""""""""""""""""""""""""""""""""""""""

#bowtie_out_statistics.log
#unfiltered_unmapped_statistics.log
#filtered_bam_statistics.log
#filtered_unmapped_statistics.log

#Pileups as bedgraphs and bigwigs : (updated FILE NAMES in the very very end of the script --> )
#""""""""""""""""""""""""""""""""""

#bowtie_out_mapped_pileup.bdg (proper pairs - no duplicates, no ploidy regions)
#bowtie_out_mapped_pileup.bw

#bowtie_out_unmapped_pileup.bdg (mapped, paired in sequencing, not mapped as proper pair - no duplicates, no ploidy regions)
#bowtie_out_unmapped_pileup.bw

#filtered_pileup.bdg (proper pairs - no duplicates, no ploidy regions)
#filtered_pileup.bw

#LEFT_pileup.bdg (footprint, left read, leftmost 1base) (proper pairs - no duplicates, no ploidy regions)
#LEFT_pileup.bw
#RIGHT_pileup.bdg (footprint, right read, rightmost 1base) (proper pairs - no duplicates, no ploidy regions)
#RIGHT_pileup.bw

#filtered_unmapped_pileup.bdg (mapped, paired in sequencing, not mapped as proper pair - no duplicates, no ploidy regions)
#filtered_unmapped_pileup.bw

#-----------------------------------------------------------------

# Runtime error handling :

# EXITS the script if bowtie_output.bam or filtered.bam (the to-be-outputted) files are not generated properly

# ECHOes "could not generate visualisation for THIS FILE" if other temporary bam files are not generated properly
# In that case never enters the visualisation generation.

# Does not check for file contents - if file exists, is readable and is not empty, it is assumed to be "good file".

#echo "Code by Jelena Telenius, 24/10/2013"


#----------------------------------------------------------------

# The code begins !

#------------------------------------------

# Printing the script name to output logs.
echo "--------------------------"
echo "$0"
which $0
echo
echo "--------------------------" >&2
echo "$0" >&2
which $0 >&2
echo >&2

# Loading subroutines in ..

echo "Loading subroutines in .."

PipeTopPath="$( which $0 | sed 's/\/bin\/mainScripts\/afterBowtieMapping.sh$//' )"

BashHelpersPath="${PipeTopPath}/bin/bashHelpers"

# PRINTING TO LOG AND ERROR FILES
. ${BashHelpersPath}/logFilePrinter.sh
# TESTING FILE EXISTENCE
. ${BashHelpersPath}/fileTesters.sh
# SETTING THE GENOME BUILD PARAMETERS
. ${BashHelpersPath}/genomeSetters.sh
# SETTING THE BLACKLIST FILTERING PARAMETERS
. ${BashHelpersPath}/blacklistSetters.sh

#------------------------------------------

# Calling in the CONFIGURATION script and its default setup :

echo "Calling in the conf/genomeBuildSetup.sh script and its default setup .."

confFolder="${PipeTopPath}/conf"

supportedGenomes=()
BOWTIE1=()
BOWTIE2=()
UCSC=()
genomesWhichHaveBlacklist=()
BLACKLIST=()

. ${confFolder}/genomeBuildSetup.sh

setGenomeLocations

echo 
echo "Supported genomes : "
for g in $( seq 0 $((${#supportedGenomes[@]}-1)) ); do echo -n "${supportedGenomes[$g]} "; done
echo 
echo

echo 
echo "Blacklist filtering available for these genomes : "
for g in $( seq 0 $((${#genomesWhichHaveBlacklist[@]}-1)) ); do echo -n "${genomesWhichHaveBlacklist[$g]} "; done
echo 
echo 

#------------------------------------------

# From where to call the main scripts operating from the mainScripts folder..

PipePath="${PipeTopPath}/bin/mainScripts"

#------------------------------------------

# From where to call the helper perl scripts..

PerlHelpersPath="${PipeTopPath}/bin/perlHelpers"

#------------------------------------------

echo
echo "PipeTopPath ${PipeTopPath}"
echo "PipePath ${PipePath}"
echo "confFolder ${confFolder}"
echo "BashHelpersPath ${BashHelpersPath}"
echo "PerlHelpersPath ${PerlHelpersPath}"
echo

#------------------------------------------


genomeName=""
MERGE_FP=1
CONTIG_FP=2
DEPTH_FP=3
MERGE_P=1
CONTIG_P=2
DEPTH_P=3
PLOIDYFILTER=1
SINGLEEND=0
PEAKCALL=0
FOOTPRINT=0
memoryMegas=265

WIN=10
WINTRACK=20
INCTRACK=2
saveUnpaired=0
saveUnfiltered=0
PipePath="/t1-home/molhaem2/telenius/Jelena_DNase_pipe_VS_100_180215"

OPTS=`getopt -o g:,d:,m:,c:,w:,u:,U:,p:,D:,M:,C:,f:,e:,P:,F:,W:,I:,b: -- "$@"`
if [ $? != 0 ]
then
    exit 1
fi

eval set -- "$OPTS"

while true ; do
    case "$1" in
        -g) genomeName=$2 ; shift 2;;
        -m) MERGE_FP=$2 ; shift 2;;
        -c) CONTIG_FP=$2 ; shift 2;;
        -d) DEPTH_FP=$2 ; shift 2;;
        -w) WIN=$2 ; shift 2;;
        -u) saveUnpaired=$2 ; shift 2;;
        -U) saveUnfiltered=$2 ; shift 2;;
        -p) PipePath=$2 ; shift 2;;
        -M) MERGE_P=$2 ; shift 2;;
        -C) CONTIG_P=$2 ; shift 2;;
        -D) DEPTH_P=$2 ; shift 2;;
        -f) PLOIDYFILTER=$2 ; shift 2;;
        -e) SINGLEEND=$2 ; shift 2;;
        -P) PEAKCALL=$2 ; shift 2;;
        -F) FOOTPRINT=$2 ; shift 2;;
        -W) WINTRACK=$2 ; shift 2;;
        -I) INCTRACK=$2 ; shift 2;;
        -b) memoryMegas=$2 ; shift 2;;
        --) shift; break;;
    esac
done

MERGE=${MERGE_FP}
CONTIG=${CONTIG_FP}
DEPTH=${DEPTH_FP}

setUCSCgenomeSizes

echo "Starting run with parameters :"
echo "genomeName ${genomeName}"
echo "WIN ${WIN}"
echo "saveUnpaired ${saveUnpaired}"
echo "saveUnfiltered ${saveUnfiltered}"
echo "PipePath ${PipePath}"
echo "MERGE_P ${MERGE_P}"
echo "CONTIG_P ${CONTIG_P}"
echo "DEPTH_P ${DEPTH_P}"
echo "BLACKLISTFILTER ${PLOIDYFILTER}"
echo "SINGLEEND ${SINGLEEND}"
echo "PEAKCALL ${PEAKCALL}"
echo "FOOTPRINT ${FOOTPRINT}"
echo

#-If-it-is-ONLY-for-peak-call-and/or-windowing-purposes..---------------------------------
if [ -d "Rerun" ]; then

    if [ -s "Rerun/LEFT_sortedTempFile.bed" ]; then

    MERGE=${MERGE_FP}
    CONTIG=${CONTIG_FP}
    DEPTH=${DEPTH_FP}

    mv -f Rerun/LEFT_sortedTempFile.bed sortedTempFile.bed
    
    ThisEnd="LEFT"
    MERGEFIRST=1
    doPeakCall
    
    WindowInput="sortedTempFile.bed"
    WindowOutputName="LEFT_window.bdg"
    WINDOW=$((${WIN}*4))
    INCREMENT=$((${WIN}*2))
    doWindowing
    rm -f sortedTempFile.bed
    
    thisBDGfile="LEFT_window.bdg"
    DoRounding=1
    doBigWigging
    
    fi
    
    if [ -s "Rerun/RIGHT_sortedTempFile.bed" ]; then
    
    mv -f Rerun/RIGHT_sortedTempFile.bed sortedTempFile.bed   
    
    ThisEnd="RIGHT"
    MERGEFIRST=1
    doPeakCall
    
    WindowInput="sortedTempFile.bed"
    WindowOutputName="RIGHT_window.bdg"
    WINDOW=$((${WIN}*4))
    INCREMENT=$((${WIN}*2))
    doWindowing
    rm -f sortedTempFile.bed 
    
    thisBDGfile="RIGHT_window.bdg"
    DoRounding=1
    doBigWigging
    
    fi
    
   
    #---------------------------------------------------------
    # GENERATING THE PEAK CALL FOR READ DEPTH
    #---------------------------------------------------------
    
    if [ -s "Rerun/READ_sortedTempFile.bed" ]; then
    
    mv -f Rerun/READ_sortedTempFile.bed sortedTempFile.bed
    
    MERGE=${MERGE_P}
    CONTIG=${CONTIG_P}
    DEPTH=${DEPTH_P}
    
    ThisEnd="READ"
    MERGEFIRST=0
    doPeakCall
    
    # Sort the peak file ..
    mv -f READ_peakCall.bdg intoSorting.txt
    sortParams='-k1,1 -k2,2n'
    sortIn1E6bunches   
    # If all went well, we delete original file. If not, we complain here (but will not die).
    sortResultInfo
    rm -f intoSorting.txt
    mv -f TEMPsortedMerged.txt READ_peakCall.bdg
    
    # And also generating the peak call quantitative file
    thisIsWhereIam=$( pwd )
    printThis="Starting to sort READ_peakCall_quant BIG time - will save temporary files in ${thisIsWhereIam}"
    printToLogFile
    # Before we use bedtools/2.x we cannot have sorted input to this. So this is possible memory peak ..
    bedtools coverage -counts -a sortedTempFile.bed -b READ_peakCall.bdg | cut -f 1,2,3,5 > intoSorting.txt
    rm -f sortedTempFile.bed
    
    sortParams='-k1,1 -k2,2n'
    sortIn1E6bunches
  # needs these to be set :
  # thisIsWhereIam=$( pwd )
  # sortParams="-k1,1 -k2,2n"  or sortParams="-n" etc
  # input in intoSorting.txt
  # outputs TEMPsortedMerged.txt
    
    # If all went well, we delete original file. If not, we complain here (but will not die).
   sortResultInfo
   rm -f intoSorting.txt
    
    mv -f TEMPsortedMerged.txt READ_peakCall_quant.bdg
    
    thisBDGfile="READ_peakCall.bdg"
    DoRounding=0
    doBigWigging
    
    thisBDGfile="READ_peakCall_quant.bdg"
    DoRounding=0
    doBigWigging
    
    echo -e "chr\tstart\tstop\treadDepth\tpeakWidth" > heading.txt
    cat READ_peakCall_quant.bdg | cut -f 2,3 | awk '{ print $2-$1 }' > TEMPwidth.txt
    paste READ_peakCall_quant.bdg TEMPwidth.txt > TEMP_intoMIG.txt
    cat heading.txt TEMP_intoMIG.txt > forMigging.txt
    rm -f TEMPwidth.txt TEMP_intoMIG.txt
    
    perl ${PerlHelpersPath}/data2gff.pl --data forMigging.txt --output PeakCall.gff
    rm -f heading.txt forMigging.txt
    
    fi
    
    rmdir Rerun

else
#-----------------------------------------------------------------------------------------

    printThis="Sam --> Bam transform, bowtie run statistics \ngenerating filtered.bam, bigwig and bedgraph files"
    printNewChapterToLogFile
    
    printThis="${genomeName} GENOME"
    printToLogFile
    
    echo "We are in folder "
    pwd 
    echo
    
if [ -e "./bowtie_out.bam" ] && [ -r "./bowtie_out.bam" ] && [ -s "./bowtie_out.bam" ] ; then  
   
    if [ "$SINGLEEND" -eq 0 ] ; then
    
    echo "Generating filtered bam, and visualisation files.."
    echo
    echo "ONLY PROPER PAIRS to all visualisation files (samtools view -b -f 2 in.bam)"
    echo "Exception the unpaired plots, which are, as mentioned UNMAPPED."
    echo
    echo "BEDGRAPH generation (bedtools bamtobed -i in.bam -bedpe | cut -f 1,2,6 | bedtools genomecov -bg -i stdin -g ${genomeName}.chrom.sizes > pileup.bdg)"
    echo "BIGWIG generation (bedGraphToBigWig pileup.bdg ${genomeName}_sizes.txt pileup.bw)"
    echo
    echo "FILTERED.BAM generation (only proper pairs --> no duplicates --> no BLACKLISTED regions):"
    echo "1) ONLY PROPER PAIRS (samtools view -b -f 2 in.bam)"
    echo "2) NO DUPLICATES (samtools rmdup in.bam out.bam)"
    echo "3) NO BLACKLISTED REGIONS (bedtools pairtobed -abam in.bam -b blacklistedRegions.bed -type neither > filtered.bam)"
    echo
    echo "UNMAPPED_FILTERED.BAM generation (only reads which failed to map in paired run --> no duplicates --> no BLACKLISTED regions):"
    echo "1) SINGLE-END MAPPED READS"
    echo "2) NO DUPLICATES (samtools rmdup in.bam out.bam)"
    echo "3) NO BLACKLISTED REGIONS (bedtools pairtobed -abam in.bam -b blacklistedRegions.bed -type neither > filtered.bam)"
    echo

    else
    
    echo "Generating visualisation files.."
    echo
    
    fi
    
    
#---------STATISTICS-for-original-file----------------------------------------------------------------------------

    samtools view -H bowtie_out.bam > original_heading.sam

    samtools sort bowtie_out.bam bowtie_out_Sorted
    samtools index bowtie_out_Sorted.bam
    testedFile="./bowtie_out_Sorted.bam"
    doTempFileTesting
    thisSortedBAMfile="bowtie_out_Sorted.bam"
    doStatisticsFile
    rm -f bowtie_out_Sorted*
    
    
#---------UNFILTERED-MAPPED-----filtering-nonmapped-or-unpaired-reads----------------------------------------------------------------------------

if [ "$SINGLEEND" -eq 0 ] ; then

    # FILE WITH PROPER PAIRS ONLY
    echo FILE WITH PROPER PAIRS ONLY >&2
    samtools view -b -f 2 bowtie_out.bam > mapped.bam
    
    testedFile="./mapped.bam"
    doTempFileTesting
    
else
# Single end gets rid of "nonmapped" reads (and gets named "paired"):
    echo FILE WITH MAPPED READS ONLY >&2
    samtools view -b -F 4 bowtie_out.bam > mapped.bam
fi


# Whether this file gets saved or not.. (flag values 10 and 11 save this file)
if [ "${saveUnfiltered}" -eq 1 ] || [ "${saveUnfiltered}" -eq 0 ]; then
rm -f bowtie_out.bam
fi

#-STATISTICS------SORTING-and-INDEXING--------------------------------------------------------------

    samtools sort mapped.bam mapped_Sorted
    samtools index mapped_Sorted.bam
    testedFile="./mapped_Sorted.bam"
    doTempFileTesting
    thisSortedBAMfile="mapped_Sorted.bam"
    doStatisticsFile
    rm -f mapped_Sorted.bam*

#---------UNFILTERED-MAPPED------filtering-unmapped-and-unpaired-out-first---------------------------------------------------------------------
    
    echo UNFILTERED bam VISUALISATION >&2
    thisIsWhereIam=$( pwd )
    printThis="Starting to sort mapped.bam BIG time - will save temporary files in ${thisIsWhereIam}"
    printToLogFile
    
    # First extracting bed if we have paired end :
    if [ "$SINGLEEND" -eq 0 ] ; then

        bedtools bamtobed -i mapped.bam -bedpe | cut -f 1,2,6 > intoSorting.txt
    
    else
    # Single end :
    
        if [ ! -r "./mapped.bam" ] || [ ! -e "./mapped.bam" ] || [ ! -f "./mapped.bam" ]; then
          echo "Temporary file not found or empty file : mapped.bam" >&2
          echo "No visualisation generated for UNFILTERED bam"
        else
        
        bedtools bamtobed -i mapped.bam | cut -f 1,2,3 > intoSorting.txt
        
        fi
    fi
    
    sortParams='-k1,1 -k2,2n'
    sortIn1E6bunches
  # needs these to be set :
  # thisIsWhereIam=$( pwd )
  # sortParams="-k1,1 -k2,2n"  or sortParams="-n" etc
  # input in intoSorting.txt
  # outputs TEMPsortedMerged.txt
    
    # If all went well, we delete original file. If not, we complain here (but will not die).
   sortResultInfo
   rm -f intoSorting.txt
   
   mv -f TEMPsortedMerged.txt forCoverages.bed
   

#---------UNFILTERED-MAPPED------if-we-DO-WINDOW-or-DO-NOT-window-over-regions---------------------------------------------------------------------

     
        bedtools genomecov -bg -i forCoverages.bed -g ${genomeBuild} > bowtie_out_mapped_pileup.bdg
        
        testedFile="bowtie_out_mapped_pileup.bdg"
        doTempFileInfo
       
    if [ "${WINTRACK}" -ne 0 ] ; then 
        
        WindowInput="forCoverages.bed"
        WindowOutputName="bowtie_out_mapped_window.bdg"
        WINDOW=${WINTRACK}
        INCREMENT=${INCTRACK}
        doWindowing
        
        testedFile="bowtie_out_mapped_window.bdg"
        doTempFileInfo
        
    fi
    
    rm -f forCoverages.bed

    #echo "--------------------------------------------------"
    #ls -lh | cut -d " " -f 1,2,3,4 --complement  >&2
    
    
#---------PLOIDY-FILTER------------------------------------------------------------------------------------

 
    # First ploidy filter (blacklisted regions filter) - carefully, as the original order should be preserved (same pairs as in the original output file ! )
    
    ploidyPath=""
    weHavePloidyFile=0
    if [ "${PLOIDYFILTER}" -eq 1 ] ; then
    setPloidyPath
    else    
    printThis="Skipping filtering BLACKLISTED regions ( blacklisted regions filtering was not requested ) .."
    printToLogFile
    fi
    
    # The above will set weHavePloidyFile=0 , if the genome wasn't listed in the genomes having blacklist file.
    
    if [ "${weHavePloidyFile}" -eq 1 ] ; then
    # Now we know, that weHavePloidyFile is only 1 , if we actually have a genome which has blacklist file.
    # Proceeding to filtering, then !
    
    printThis="Filtering BLACKLISTED regions .."
    printToLogFile

        if [ "$SINGLEEND" -eq 0 ] ; then
        bedtools pairtobed -abam mapped.bam -b "${ploidyPath}" -type neither > noPloidyMapped.bam
        else  
        bedtools intersect -v -abam mapped.bam -b "${ploidyPath}" > noPloidyMapped.bam
        fi  
    
        # Whether this file gets saved or not.. (flag values 01 and 11 save this file)
        if [ "${saveUnfiltered}" -eq 0 ] || [ "${saveUnfiltered}" -eq 10 ]; then
        rm -f mapped.bam
        fi
    
    else
    
    printThis="Skipping filtering BLACKLISTED regions ( blacklisted regions file not available for genome ${genomeName} ) .."
    printToLogFile
    
        if [ "${saveUnfiltered}" -eq 0 ] || [ "${saveUnfiltered}" -eq 10 ]; then
        mv -f mapped.bam noPloidyMapped.bam
        else
        cp mapped.bam noPloidyMapped.bam
        fi
    
    fi

    testedFile="./noPloidyMapped.bam"
    doTempFileTesting
    
#---------DUPLICATE-FILTER------------------------------------------------------------------------------------

    samtools sort noPloidyMapped.bam noPloidyMapped_Sorted
    
    testedFile="./noPloidyMapped_Sorted.bam"
    doTempFileTesting
    rm -f noPloidyMapped.bam

    # FILTERED BAM GENERATION
    echo FILTERED BAM GENERATION >&2
    if [ "$SINGLEEND" -eq 0 ] ; then
    samtools rmdup noPloidyMapped_Sorted.bam filtered_Sorted.bam
    else
    # Single end file gets named "paired"
    samtools rmdup -s noPloidyMapped_Sorted.bam filtered_Sorted.bam
    fi    
    
    testedFile="./filtered_Sorted.bam"
    doTempFileTesting
    
    rm -f noPloidyMapped_Sorted.bam
    
#-----------------------------------------------------------
    
    
    if [ ! -r "./filtered_Sorted.bam" ] || [ ! -e "./filtered_Sorted.bam" ] || [ ! -f "./filtered_Sorted.bam" ]; then
      echo "Temporary file not found or empty file : filtered_Sorted.bam" >&2
      echo "No visualisations or footprints generated FILTERED bam"
    else    

    # FILTERED BAM STATISTICS
    
    # Already sorted, no need to sort again..
    #samtools sort filtered.bam filtered_Sorted
    samtools index filtered_Sorted.bam
    
    thisSortedBAMfile="filtered_Sorted.bam"
    doStatisticsFile

    # FILTERED bam VISUALISATION
    
    echo FILTERED bam VISUALISATION >&2
   
    
   #--ORIGINAL-ORDERING--------------------------------------------------------
   
   samtools view filtered_Sorted.bam > filtered_SortedNoheading.sam
   # using the previously generated original_heading.sam (which says "unsorted" in the sorting field)
   mv -f original_heading.sam TEMP_heading.sam
   #samtools view -H filtered_Sorted.bam > TEMP_heading.sam
   
    testedFile="./filtered_SortedNoheading.sam"
    doTempFileTesting
    
    # Move the bams to be TEMPORARY STORED - only saved if we fail in the sorting to original (so that user can salvage the crashed run outside-pipeline )
    mv -f filtered_Sorted.bam filtered_coordSorted_SAFETYCOPY.bam
    mv -f filtered_Sorted.bam.bai filtered_coordSorted_SAFETYCOPY.bam.bai
   
   # Fixing over-greedy sort 09Sep2016 (working diary 39).
   # resulting to add these :
   # -S 7G  (to avoid bug in "only compressing the first file" - described here :
   #          http://unix.stackexchange.com/questions/275501/gnu-sort-compress-program-compressing-only-first-temporary
   # (better compression than gzip)
   # --temporary-directory=${rundir} or shorter -T .
   
   thisIsWhereIam=$( pwd )
   printThis="Starting to sort filtered_SortedNoheading.sam BIG time - will save temporary files in ${thisIsWhereIam}"
   printToLogFile
   cat filtered_SortedNoheading.sam | sed 's/^.*CO:Z://' | sed 's/\s*//' | paste - filtered_SortedNoheading.sam > intoSorting.txt

   # If all went well, we delete original file.
    countBefore=$(($( cat filtered_SortedNoheading.sam | grep -c "" )))
    countAfter=$(($( cat intoSorting.txt | grep -c "" )))
    
    if [ "${countBefore}" -ne "${countAfter}" ]
    then
    
    echo "Preparing filtered_SortedNoheading.sam for sorting to original order FAILED. " >&2
    echo "Original file had ${countBefore} data lines - prepared file had only ${countAfter} lines." >&2
    echo "EXITING!!" >&2
    exit 1

    fi   
   
    rm -f filtered_SortedNoheading.sam
  
  
  # Now, sorting.
  
  sortParams='-n'
  sortIn1E6bunches
  # needs these to be set :
  # thisIsWhereIam=$( pwd )
  # sortParams="-k1,1 -k2,2n"  or sortParams="-n" etc
  # input in intoSorting.txt
  # outputs TEMPsortedMerged.txt

   # If all went well, we delete original file. If not, we die here.
   sortResultTester
   rm -f intoSorting.txt
    
    # We combine to final filtered file (now back in original order)
    
    cut -f 1 --complement TEMPsortedMerged.txt | cat TEMP_heading.sam -  > filtered.sam
   
    testedFile="./filtered.sam"
    doTempFileTesting
    
    rm -f TEMP_heading.sam TEMPsortedMerged.txt 
   
   samtools view -bS filtered.sam > filtered.bam
   
    testedFile="./filtered.bam"
    doTempFileTesting
    rm -f filtered.sam
    rm -f filtered_coordSorted_SAFETYCOPY.bam*
    
    thisIsWhereIam=$( pwd )
    printThis="Starting to sort filtered.bam BIG time - will save temporary files in ${thisIsWhereIam}"
    printToLogFile
    
    if [ "$SINGLEEND" -eq 0 ] ; then
    bedtools bamtobed -i filtered.bam -bedpe | cut -f 1,2,6 > intoSorting.txt 
    else
    bedtools bamtobed -i filtered.bam | cut -f 1,2,3 > intoSorting.txt
    fi
    
    sortParams='-k1,1 -k2,2n'
    sortIn1E6bunches
  # needs these to be set :
  # thisIsWhereIam=$( pwd )
  # sortParams="-k1,1 -k2,2n"  or sortParams="-n" etc
  # input in intoSorting.txt
  # outputs TEMPsortedMerged.txt
    
    # If all went well, we delete original file. If not, we complain here (but will not die).
   sortResultInfo
   rm -f intoSorting.txt
   
   mv -f TEMPsortedMerged.txt filtered_sorted.bed
    
    
    bedtools genomecov -bg -i filtered_sorted.bed -g ${genomeBuild} > filtered_pileup.bdg
    
    testedFile="filtered_pileup.bdg"
    doTempFileInfo
    
    if [ "${WINTRACK}" -ne 0 ]; then 
    
    WindowInput="filtered_sorted.bed"
    WindowOutputName="filtered_window.bdg"
    WINDOW=${WINTRACK}
    INCREMENT=${INCTRACK}
    doWindowing
    
    testedFile="filtered_window.bdg"
    doTempFileInfo
   
    fi
    
    #---------------------------------------------------------
    # GENERATING THE PEAK CALL 
    #---------------------------------------------------------
    
    mkdir ForPeakCallAndFootPrintRerun
    
    mv filtered_sorted.bed sortedTempFile.bed
    
    if [ "${PEAKCALL}" -eq 1 ] ; then
    
    MERGE=${MERGE_P}
    CONTIG=${CONTIG_P}
    DEPTH=${DEPTH_P}
    
    ThisEnd="READ"
    MERGEFIRST=0
    doPeakCall
    
    # And also generating the peak call quantitative file
    thisIsWhereIam=$( pwd )
    printThis="Starting to sort READ_peakCall_quant file BIG time - will save temporary files in ${thisIsWhereIam}"
    printToLogFile
    bedtools coverage -counts -a sortedTempFile.bed -b READ_peakCall.bdg | cut -f 1,2,3,5 > intoSorting.txt
    
    sortParams='-k1,1 -k2,2n'
    sortIn1E6bunches
  # needs these to be set :
  # thisIsWhereIam=$( pwd )
  # sortParams="-k1,1 -k2,2n"  or sortParams="-n" etc
  # input in intoSorting.txt
  # outputs TEMPsortedMerged.txt
    
    # If all went well, we delete original file. If not, we complain here (but will not die).
   sortResultInfo
   rm -f intoSorting.txt
   
   mv -f TEMPsortedMerged.txt READ_peakCall_quant.bdg
    
    
    # end of this if : if [ "${PEAKCALL}" -eq 1 ] ; then
    fi
    
    mv -f sortedTempFile.bed ForPeakCallAndFootPrintRerun/READ_sortedTempFile.bed
    #tree
    
    #---------------------------------------------------------
    # GENERATING THE FOOTPRINTS
    #---------------------------------------------------------
    
    # FILTERED bam FOOTPRINT
    echo FILTERED bam FOOTPRINT >&2
    
    # We produce pileups for i) LEFTMOST end of left reads (1bp), and ii) RIGHTMOST end of right reads (1bp) - to monitor WHERE EXACTLY our DNaseI cuts were !
    
    # Here we end up fixing the RIGHTMOST end of right reads (which we of course want to pile) 1base to left, because :
    #http://samtools.sourceforge.net/SAMv1.pdf
    
    #1-based coordinate system 

    #A coordinate system where the first base of a sequence is one. 
    #In this coordinate system, a region is specified by a closed interval. 
    #For example, the region between the 3rd and the 7th bases inclusive is [3;7]. 
    #The SAM, VCF, GFF and Wiggle formats are using the 1-based coordinate system.

    #0-based coordinate system 

    #A coordinate system where the first base of a sequence is zero. 
    #In this coordinate system, a region is specified by a half-closed-half-open interval. 
    #For example, the region between the 3rd and the 7th bases inclusive is [2;7). 
    #The BAM, BCFv2, BED, and PSL formats are using the 0-based coordinate system.
        
    # Left reads :
    echo Left reads >&2
    
    # Paired end first..
    if [ "$SINGLEEND" -eq 0 ] ; then
    #Awk, because the coordinate is start position, and if end position is the same base, then it needs to be +1, as the region is [start,stop[ -i.e. half open
    bedtools bamtobed -i filtered.bam -bedpe | cut -f 1,2 > tempFile.txt
    # Single end.. 
    else

    # Normal single end
    samtools view -b filtered.bam | bedtools bamtobed -i stdin | cut -f 1,2,3,6 > tempFileCombo.txt
    
    # which map to (+) strand
    cat tempFileCombo.txt | awk '$4 == "+" {print $0}' > tempLeft.txt
    # which map to (-) strand
    cat tempFileCombo.txt | awk '$4 == "-" {print $0}' | cut -f 1,3 > tempRight.txt
    
    rm -f tempFileCombo.txt
    
    mv -f tempLeft.txt tempFile.txt
    
    echo "Generated temp files :"
    ls -lht
    
    fi
    
    # Both single and paired end continue the same way..
    
    testedFile="tempFile.txt"
    doTempFileInfo
    
    cut -f 2 tempFile.txt | awk '{s=$1+1} {print s}' > tempCol3.txt
    thisIsWhereIam=$( pwd )
    printThis="Starting to sort tempFile.txt BIG time - will save temporary files in ${thisIsWhereIam}"
    printToLogFile
    cut -f 1,2 tempFile.txt | paste - tempCol3.txt > intoSorting.txt
    rm -f tempFile.txt tempCol3.txt
    
    sortParams='-k1,1 -k2,2n'
    sortIn1E6bunches
  # needs these to be set :
  # thisIsWhereIam=$( pwd )
  # sortParams="-k1,1 -k2,2n"  or sortParams="-n" etc
  # input in intoSorting.txt
  # outputs TEMPsortedMerged.txt
    
    # If all went well, we delete original file. If not, we complain here (but will not die).
   sortResultInfo
   rm -f intoSorting.txt
    
   mv -f  TEMPsortedMerged.txt sortedTempFile.bed
    
    
    #--------------------------------------------
    # SETTING WINDOW SIZE AND INCREMENT HERE
    #--------------------------------------------
    
    if [ "${FOOTPRINT}" -eq 1 ] ; then
      
    echo
    echo "head LEFT_pileup.bdg"
    head LEFT_pileup.bdg
    echo
    
    WindowInput="sortedTempFile.bed"
    WindowOutputName="LEFT_window.bdg"
    WINDOW=$((${WIN}*4))
    INCREMENT=$((${WIN}*2))
    doWindowing
    
    # End of this if : if [ "${FOOTPRINT}" -eq 1 ] ; then
    fi
    
    mv -f sortedTempFile.bed ForPeakCallAndFootPrintRerun/LEFT_sortedTempFile.bed
   
    # Right reads (including the awk oneliner moving the position 1base to left):
    echo Right reads >&2
    
    # Paired end first..
    if [ "$SINGLEEND" -eq 0 ] ; then
    bedtools bamtobed -i filtered.bam -bedpe | cut -f 1,6 > tempFile.txt
    # Single end..
    else

    mv -f tempRight.txt tempFile.txt
    
    fi
    
    # Both single and paired end continue the same way..
    
    testedFile="tempFile.txt"
    doTempFileInfo
    
    cut -f 1 tempFile.txt > tempCol1.txt
    cut -f 2 tempFile.txt | awk '{s=$1-1} {print s}' > tempCol2.txt
    cut -f 2 tempFile.txt > tempCol3.txt
    thisIsWhereIam=$( pwd )
    printThis="Starting to sort sortedTempFile BIG time - will save temporary files in ${thisIsWhereIam}"
    printToLogFile
    
    testedFile="tempCol1.txt"
    doTempFileInfo
    testedFile="tempCol2.txt"
    doTempFileInfo
    testedFile="tempCol3.txt"
    doTempFileInfo
    
    paste tempCol1.txt tempCol2.txt tempCol3.txt > intoSorting.txt
    
    sortParams='-k1,1 -k2,2n'
    sortIn1E6bunches
  # needs these to be set :
  # thisIsWhereIam=$( pwd )
  # sortParams="-k1,1 -k2,2n"  or sortParams="-n" etc
  # input in intoSorting.txt
  # outputs TEMPsortedMerged.txt
    
    # If all went well, we delete original file. If not, we complain here (but will not die).
   sortResultInfo
   rm -f intoSorting.txt
   
   mv -f TEMPsortedMerged.txt sortedTempFile.bed
    

    testedFile="sortedTempFile.bed"
    doTempFileInfo
    
    rm -f tempFile.txt tempCol1.txt tempCol2.txt tempCol3.txt
    
    if [ "${FOOTPRINT}" -eq 1 ] ; then
    
    echo
    echo "head RIGHT_pileup.bdg"
    head RIGHT_pileup.bdg
    echo
    
    WindowInput="sortedTempFile.bed"
    WindowOutputName="RIGHT_window.bdg"
    WINDOW=$((${WIN}*4))
    INCREMENT=$((${WIN}*2))
    doWindowing
    
    # End of this if : if [ "${FOOTPRINT}" -eq 1 ] ; then
    fi
    
    mv -f sortedTempFile.bed ForPeakCallAndFootPrintRerun/RIGHT_sortedTempFile.bed
    #tree
    
    #---------------------------------------------------------
    # Packing the temporary files..
    
    cd ForPeakCallAndFootPrintRerun
    thisBEDfile="LEFT_sortedTempFile.bed"
    doBigBedding
    rm -f LEFT_sortedTempFile.bed
    thisBEDfile="RIGHT_sortedTempFile.bed"
    doBigBedding
    rm -f RIGHT_sortedTempFile.bed
    thisBEDfile="READ_sortedTempFile.bed"
    doBigBedding
    rm -f READ_sortedTempFile.bed
    cd ..

    #---------------------------------------------------------
    
    #Exiting visualisation 'if filtered.bam file exists'    
    fi
    
    echo "--------------------------------------------------"

    
#---------unpaired-UNFILTERED------------------------------------------------------------------------------------

    
    # UNFILTERED unpaired VISUALISATION
    
    if [ "$SINGLEEND" -eq 0 ] ; then
    
    echo "UNFILTERED unpaired VISUALISATION - (causes samtools errors if no unpaired reads found)" >&2
    
    if [ "$(($( samtools view -c singleEnd_bowtie_READ1.bam )))" -eq 0 ]; then
        echo "No only-single-end-mapped reads in READ1 ! "
    else
        
        
        
        
        samtools sort singleEnd_bowtie_READ1.bam singleEnd_READ1_Sorted
        testedFile="./singleEnd_READ1_Sorted.bam"
        doTempFileTesting
        
        # If we want to save this stage as bam file (save when flag value 10 or 11)
        if [ "${saveUnpaired}" -eq 1 ] || [ "${saveUnpaired}" -eq 0 ]; then
            rm -f singleEnd_bowtie_READ1.bam
        fi
        
        samtools index singleEnd_READ1_Sorted.bam        
        unpairedSorted="singleEnd_READ1_Sorted.bam"
        unpairedBdg  
    fi
    
    
    if [ "$(($( samtools view -c singleEnd_bowtie_READ2.bam )))" -eq 0 ]; then
        echo "No only-single-end-mapped reads in READ2 ! "
    else
        samtools sort singleEnd_bowtie_READ2.bam singleEnd_READ2_Sorted
        testedFile="./singleEnd_READ2_Sorted.bam"
        doTempFileTesting
        
        # If we want to save this stage as bam file (save when flag value 10 or 11)
        if [ "${saveUnpaired}" -eq 1 ] || [ "${saveUnpaired}" -eq 0 ]; then
            rm -f singleEnd_bowtie_READ2.bam
        fi
        
        
        samtools index singleEnd_READ2_Sorted.bam
        unpairedSorted="singleEnd_READ2_Sorted.bam"
        unpairedBdg  
    fi
    
    fi
    
    
#---------BIGWIG-GENERATION------------------------------------------------------------------------------------

    doBigWigGeneration

#-------------------------------------------    
# If the file bowtie_out.sam was not found :
else
    echo "file bowtie_out.bam was not found - aborting bowtie statistics, bam filtering and bigwig generation" >&2
fi

fi