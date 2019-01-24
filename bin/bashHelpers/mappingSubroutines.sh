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

runBowtie(){
    # NEEDS THIS TO BE SET BEFORE CALL :
    # genomeName=""  - being one of these "mm9", "mm10", "hg18", "hg_19"  
    
    
    # Testing that fastq files exist !

if [ -r "./READ1.fastq" ] ; then
if [ -r "./READ2.fastq" ] || [ "${SINGLE_END}" -eq 1 ] ; then

setMparameter
setBowtieReadList

# Generating the GENOME folder
echo ""
mkdir ${genomeName}
cd ${genomeName}
printThis="Generated run directory for the genome :"
printToLogFile
pwd 
pwd  >&2

# Copying the fastq's (now as in VS14 trimming is dependent on genome build)

# If we only have ONE genome, we can MOVE the fastq files..
if [ "${#GENOMEARRAY[@]}" -eq 1 ] ; then
    mv -f ../READ1.fastq .
    if [ "${SINGLE_END}" -eq 0 ] ; then
    mv -f ../READ2.fastq .
    fi
else  
    cp ../READ1.fastq .
    if [ "${SINGLE_END}" -eq 0 ] ; then
    cp ../READ2.fastq .
    fi
fi

# Standard bowtie - to generate the "all reads" bam for reference (saving space) purposes. Only do this if input is given as FASTQ files.

if [ "${UNTRIMMED}" -ge 1 ] ; then
    runRawBowtie
fi


# Here the subroutine, which checks for all TRIMMING and FLASHING needs for the data
if [ "${TRIM}" -eq 1 ] || [ "${FLASH}" -eq 1 ]; then
    
    echoTrimInfo
    trimmingFlashing
    
    # After-trim  QC (assumes excisting READ1.fastq, READ2.fastq)
    echo "Running FastQC for trimmed (and/or) flashed reads.."
    printThis="${PipePath}/QC_and_Trimming.sh --fastqc --single ${SINGLE_END}"
    printToLogFile
    ${PipePath}/QC_and_Trimming.sh --fastqc --single ${SINGLE_END}
    
    # Changing names of fastqc folders to be "TRIMMED"
    mkdir READ1_fastqc_TRIMMED
    mv -f READ1_fastqc.html READ1_fastqc_TRIMMED/fastqc_report.html
    
    if [ "$SINGLE_END" -eq 0 ] ; then
    mkdir READ2_fastqc_TRIMMED
    mv -f READ2_fastqc.html READ2_fastqc_TRIMMED/fastqc_report.html
    fi
    
    mv -f "READ1_fastqc.zip" "READ1_fastqc_TRIMMED.zip"
    
    if [ "$SINGLE_END" -eq 0 ] ; then
    mv -f "READ2_fastqc.zip" "READ2_fastqc_TRIMMED.zip"
fi

fi
# It reads the user-set flags (trim or no trim, flash or no flash) - and trims the data accordingly, and provides fastqc as a "proof" of trimming.
# Replaces the READ1.fastq and READ2.fastq with the trimmed versions, for the subsequent part of the "runBowtie" routine below.

if [ "$SINGLE_END" -eq 0 ] ; then

printThis="Bowtie in ${genomeName} build - single end for initially unmapped reads"
printNewChapterToLogFile
    
    echo "Beginning bowtie run (first mapping paired - saving only TRULY unmapped (not excluded due -m parameter)) .."    
 
    echo  >&2
    echo "Preliminary bowtie to separate unpaired reads" >&2
    
    if [ "${BOWTIE}" -eq 2 ] ; then
    echo "bowtie2 -p ${BOWTIE_PROCESSORS} ${otherBowtie2Parameters} ${bowtieQuals} --no-discordant --no-mixed --maxins ${MAXINS} -x ${bowtieGenomeBuild} ${bowtieReadList} --un-conc UNMAPPED.fastq > /dev/null"
    bowtie2 -p ${BOWTIE_PROCESSORS} ${otherBowtie2Parameters} ${bowtieQuals} --no-discordant --no-mixed --maxins ${MAXINS} -x ${bowtieGenomeBuild} ${bowtieReadList} --un-conc UNMAPPED.fastq > /dev/null
    else
    echo "bowtie -p ${BOWTIE_PROCESSORS} ${otherBowtie1Parameters} --chunkmb "${BOWTIEMEMORY}" ${bowtieQuals} ${mParameter} --best --strata --maxins ${MAXINS} --sam ${bowtieGenomeBuild} ${bowtieReadList} --un UNMAPPED.fastq --max M_FILTERED.fastq > /dev/null"
    bowtie -p ${BOWTIE_PROCESSORS} ${otherBowtie1Parameters} --chunkmb "${BOWTIEMEMORY}" ${bowtieQuals} ${mParameter} --best --strata --maxins ${MAXINS} --sam ${bowtieGenomeBuild} ${bowtieReadList} --un UNMAPPED.fastq --max M_FILTERED.fastq > /dev/null
    fi
    
    if [ "${BOWTIE}" -eq 2 ] ; then
       mv -f UNMAPPED.1.fastq UNMAPPED_1.fastq
       mv -f UNMAPPED.2.fastq UNMAPPED_2.fastq
    fi
    
    rm -f M_FILTERED*.fastq
    
    if [ -s "UNMAPPED_1.fastq" ] ; then

    echo  >&2
    echo "Beginning bowtie run (single end run for the unmapped READ1) - outputting run command after completion.."
    echo "Beginning bowtie run (single end run for the unmapped READ1)"  >&2

    if [ "${BOWTIE}" -eq 2 ] ; then
    bowtie2 -p ${BOWTIE_PROCESSORS} ${otherBowtie2Parameters} ${bowtieQuals} --maxins ${MAXINS} -x ${bowtieGenomeBuild} UNMAPPED_1.fastq  > singleEnd_bowtie_READ1.sam
    echo "bowtie2 -p ${BOWTIE_PROCESSORS} ${otherBowtie2Parameters} ${bowtieQuals} --maxins ${MAXINS} -x ${bowtieGenomeBuild} UNMAPPED_1.fastq"
    else
    bowtie -p ${BOWTIE_PROCESSORS} ${otherBowtie1Parameters} --chunkmb "${BOWTIEMEMORY}" ${bowtieQuals} ${mParameter} --sam ${bowtieGenomeBuild} UNMAPPED_1.fastq  > singleEnd_bowtie_READ1.sam
    fi
    
    if [ "${SAVEUNMAPPED}" -eq "0" ]; then
        rm -f UNMAPPED_1.fastq
    fi
    samtools view -SH singleEnd_bowtie_READ1.sam | grep bowtie
    #Making bam, printing only mapped reads
    #samtools view -F 4 -Sb -o singleEnd_bowtie_READ1.bam singleEnd_bowtie_READ1.sam
    
    samtools view -Sb -o singleEnd_bowtie_READ1.bam singleEnd_bowtie_READ1.sam

    
    testedFile="singleEnd_bowtie_READ1.bam"
    doTempFileTesting
    rm -f singleEnd_bowtie_READ1.sam
    
    fi

    if [ -s "UNMAPPED_2.fastq" ] ; then
    
    echo  >&2
    echo "Beginning bowtie run (single end run for the unmapped READ2) - outputting run command after completion.."
    echo "Beginning bowtie run (single end run for the unmapped READ2)"  >&2

    if [ "${BOWTIE}" -eq 2 ] ; then
    bowtie2 -p ${BOWTIE_PROCESSORS} ${otherBowtie2Parameters} ${bowtieQuals} --maxins ${MAXINS} -x ${bowtieGenomeBuild} UNMAPPED_2.fastq  > singleEnd_bowtie_READ2.sam
    echo "bowtie2 -p ${BOWTIE_PROCESSORS} ${otherBowtie2Parameters} ${bowtieQuals} --maxins ${MAXINS} -x ${bowtieGenomeBuild} UNMAPPED_2.fastq"
    else
    bowtie -p ${BOWTIE_PROCESSORS} ${otherBowtie1Parameters} --chunkmb "${BOWTIEMEMORY}" ${bowtieQuals} ${mParameter} --sam ${bowtieGenomeBuild} UNMAPPED_2.fastq  > singleEnd_bowtie_READ2.sam
    fi

    if [ "${SAVEUNMAPPED}" -eq "0" ]; then
        rm -f UNMAPPED_2.fastq
    fi 
    samtools view -SH singleEnd_bowtie_READ2.sam | grep bowtie
    #Making bam, printing only mapped reads
    #samtools view -F 4 -Sb -o singleEnd_bowtie_READ2.bam singleEnd_bowtie_READ2.sam
    
    samtools view -Sb -o singleEnd_bowtie_READ2.bam singleEnd_bowtie_READ2.sam

    
    testedFile="singleEnd_bowtie_READ2.bam"
    doTempFileTesting
    rm -f singleEnd_bowtie_READ2.sam
    
    fi

fi

if [ "$SINGLE_END" -eq 0 ] ; then
printThis="Bowtie in ${genomeName} build - 'ordinary bowtie'\n- i.e.reporting ONLY reads mapping in proper pairs as 'mapped' in the output bam file"
printNewChapterToLogFile
else
printThis="Bowtie in ${genomeName} build - 'ordinary bowtie'"
printNewChapterToLogFile    
fi

    echo ""
    echo "Beginning bowtie run (outputting run command after completion) .."    
    pwd 
    pwd  >&2
   
    if [ "${BOWTIE}" -eq 2 ] ; then
    bowtie2 -p ${BOWTIE_PROCESSORS} ${otherBowtie2Parameters} ${bowtieQuals} --no-discordant --no-mixed --maxins ${MAXINS} -x ${bowtieGenomeBuild} ${bowtieReadList} > bowtie_out.sam
    echo "bowtie2 -p ${BOWTIE_PROCESSORS} ${otherBowtie2Parameters} ${bowtieQuals} --no-discordant --no-mixed --maxins ${MAXINS} -x ${bowtieGenomeBuild} ${bowtieReadList}"
    else
    bowtie -p ${BOWTIE_PROCESSORS} ${otherBowtie1Parameters} --chunkmb "${BOWTIEMEMORY}" ${bowtieQuals} ${mParameter} --best --strata --maxins ${MAXINS} --sam ${bowtieGenomeBuild} ${bowtieReadList}  > bowtie_out.sam
    fi
    
    samtools view -SH bowtie_out.sam | grep bowtie
    
    generateCommentFieldAndSaveInBam

    deleteFastqLocal

    cd ..

# Ending "if files exist" test clause
else
    echo "Can not find file READ2.fastq - aborting bowtie"
fi
else
    echo "Can not find file READ1.fastq - aborting bowtie"  
fi
}

runRawBowtie(){
    
    # NEEDS THIS TO BE SET BEFORE CALL :
    # genomeName=""  - being one of these "mm9", "mm10", "hg18", "hg_19"  
    
    
    # Testing that fastq files exist !

if [ -r "./READ1.fastq" ] ; then
if [ -r "./READ2.fastq" ] || [ "${SINGLE_END}" -eq 1 ] ; then
    
setBOWTIEgenomeSizes
setMparameter

# Generating the GENOME folder
echo ""
#mkdir "ALL_READS_bam_${genomeName}"
#cd "ALL_READS_bam_${genomeName}"
pwd 
pwd  >&2

# Standard bowtie - to generate the "all reads" bam for reference (saving space) purposes.

printThis="Generating 'bam copy' of the input FASTQ data in ${genomeName} build"
printNewChapterToLogFile

    echo ""
    echo "Beginning bowtie run (outputting run command after completion) .."    
    pwd 
    pwd  >&2

    if [ "${BOWTIE}" -eq 2 ] ; then
    bowtie2 -p ${BOWTIE_PROCESSORS} ${otherBowtie2Parameters} ${bowtieQuals} --no-discordant --no-mixed --maxins ${MAXINS} -x ${bowtieGenomeBuild} ${bowtieReadList} > bowtie_out.sam
    echo "bowtie2 -p ${BOWTIE_PROCESSORS} ${otherBowtie2Parameters} ${bowtieQuals} --no-discordant --no-mixed --maxins ${MAXINS} -x ${bowtieGenomeBuild} ${bowtieReadList}"
    else
    bowtie -p ${BOWTIE_PROCESSORS} ${otherBowtie1Parameters} --chunkmb "${BOWTIEMEMORY}" ${bowtieQuals} ${mParameter} --best --strata --maxins ${MAXINS} --sam ${bowtieGenomeBuild} ${bowtieReadList} > bowtie_out.sam
    fi
    
    samtools view -SH bowtie_out.sam | grep bowtie
    
    samtools view -Sb -o UNTRIMMED_bowtie_out.bam bowtie_out.sam
    
    testedFile="UNTRIMMED_bowtie_out.bam"
    doTempFileTesting
    rm -f bowtie_out.sam
    
    samtools sort UNTRIMMED_bowtie_out.bam UNTRIMMED_bowtie_out_Sorted
    
    testedFile="UNTRIMMED_bowtie_out_Sorted.bam"
    doTempFileTesting
    #rm -f bowtie_out.bam
    
    samtools index UNTRIMMED_bowtie_out_Sorted.bam
    
    thisSortedBAMfile="UNTRIMMED_bowtie_out_Sorted.bam"
    doStatisticsFile
    
    rm -f UNTRIMMED_bowtie_out_Sorted.bam*
    
    #---------UNFILTERED-PAIRED-----filtering-nonmapped-or-unpaired-reads----------------------------------------------------------------------------
    
    if [ "$UNTRIMMED" -eq 1 ] || [ "$UNTRIMMED" -eq 11 ] ; then
        
    if [ "$SINGLE_END" -eq 0 ] ; then

    # FILE WITH PROPER PAIRS ONLY
    echo FILE WITH PROPER PAIRS ONLY >&2
    samtools view -b -f 2 UNTRIMMED_bowtie_out.bam > UNTRIMMED_onlyMapped_bowtie_out.bam
    
    # UNFILTERED bam VISUALISATION

    else
    # Single end gets rid of "nonmapped" reads (and gets named "paired"):
    echo FILE WITH MAPPED READS ONLY >&2
    samtools view -b -F 4 UNTRIMMED_bowtie_out.bam > UNTRIMMED_onlyMapped_bowtie_out.bam
    fi
    
    testedFile="UNTRIMMED_onlyMapped_bowtie_out.bam"
    doTempFileTesting
    
    samtools sort UNTRIMMED_onlyMapped_bowtie_out.bam UNTRIMMED_onlyMapped_bowtie_out_Sorted
    
    testedFile="UNTRIMMED_onlyMapped_bowtie_out_Sorted.bam"
    doTempFileTesting
    #rm -f bowtie_out.bam
    
    samtools index UNTRIMMED_onlyMapped_bowtie_out_Sorted.bam
    
    thisSortedBAMfile="UNTRIMMED_onlyMapped_bowtie_out_Sorted.bam"
    doStatisticsFile
    
    rm -f UNTRIMMED_onlyMapped_bowtie_out_Sorted.bam*

    fi
    
    # If we want to save only the filtered file :
    
    if [ "$UNTRIMMED" -eq 1  ] ; then
        
    rm -f UNTRIMMED_bowtie_out.bam
        
    fi

#    cd ..

# Ending "if files exist" test clause
else
    echo "Can not find file READ2.fastq - aborting bowtie"
fi
else
    echo "Can not find file READ1.fastq - aborting bowtie"  
fi
    
}

generateCommentFieldAndSaveInBam(){
    # Generate CO comment field in bam - this will hold the ORIGINAL READ ORDER which bowtie outputted (keeps read pairs with each others, if sorted along this flag)
    cat bowtie_out.sam | grep "^@" > bowtie_out_heading.sam
    cat bowtie_out.sam | grep -v "^@" > bowtie_out_noheading.sam
    rm -f bowtie_out.sam
    
    cut -f 1 bowtie_out_noheading.sam | cat -n | sed 's/^\s\s*//' | sed 's/^/CO:Z:/' | cut -f 1 | paste bowtie_out_noheading.sam - | cat bowtie_out_heading.sam -  > bowtie_out.sam
  
    testedFile="bowtie_out.sam"
    doTempFileTesting
    rm -f bowtie_out_noheading.sam bowtie_out_heading.sam
    
    samtools view -Sb -o bowtie_out.bam bowtie_out.sam
    
    testedFile="bowtie_out.bam"
    doTempFileTesting
    rm -f bowtie_out.sam
}

