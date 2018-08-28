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


printRunStartArrays(){
    
    echo
    echo "Genomes the data is to be mapped in this run :"
    for g in $( seq 0 $((${#GENOMEARRAY[@]} - 1)) ); do
    
        echo "${GENOMEARRAY[$g]}"

    done
    echo
    
    echo
    echo "-------------------------------------"
    echo
    echo "Ready to run ! - here printout of main for loop parameters : "
    echo
   
    for k in $( seq 0 $((${#folderList[@]} - 1)) ); do
        echo "folderList[$k]  ${folderList[$k]}"
    done    

    echo    

    for k in $( seq 0 $((${#fileList[@]} - 1)) ); do
        echo "fileList[$k]  ${fileList[$k]}"
    done    

    echo    

    for k in $( seq 0 $((${#GENOMEARRAY[@]} - 1)) ); do
        echo "GENOMEARRAY[$k]  ${GENOMEARRAY[$k]}"
    done    

    echo    
    echo "-------------------------------------"
    echo
    
}

printRunStartArraysFastq(){
    
    echo
    echo "Genomes the data is to be mapped in this run :"
    for g in $( seq 0 $((${#GENOMEARRAY[@]} - 1)) ); do
    
        echo "${GENOMEARRAY[$g]}"

    done
    echo
    
    echo
    echo "-------------------------------------"
    echo
    echo "Ready to run ! - here printout of main for loop parameters : "
    echo
   
    for k in $( seq 0 $((${#folderList[@]} - 1)) ); do
        echo "folderList[$k]  ${folderList[$k]}"
    done    

    echo    

    for k in $( seq 0 $((${#fileList1[@]} - 1)) ); do
        echo "fileList1[$k]  ${fileList1[$k]}"
    done    
    echo    

    if [ "${SINGLE_END}" -eq 0 ] ; then
      
    for k in $( seq 0 $((${#fileList2[@]} - 1)) ); do
        echo "fileList2[$k]  ${fileList2[$k]}"
    done    
    echo          
        
    fi
    
    for k in $( seq 0 $((${#GENOMEARRAY[@]} - 1)) ); do
        echo "GENOMEARRAY[$k]  ${GENOMEARRAY[$k]}"
    done    

    echo    
    echo "-------------------------------------"
    echo
    
}


runBowtieSubscript(){
    
#---------Running-bam-file-generation-and-BOWTIE-run-statistics---------------------

    echo ""
    cd ${genomeName}
    
    printThis="${PipePath}/afterBowtieMapping.sh -g ${genomeName} -m ${MERGE_FP} -c ${CONTIG_FP} -d ${DEPTH_FP} -w ${WINDOW} -u ${UNPAIR} -U ${UNFILTER} -p ${PipePath} -M ${MERGE_P} -C ${CONTIG_P} -D ${DEPTH_P} -f ${PLOIDYFILTER} -e ${SINGLE_END} -P ${PEAKCALL} -F ${FOOTPRINT} -W ${WINDOWTRACK} -I ${WINDOWINCR} -b ${BOWTIEMEMORY} "
    printToLogFile
    ${PipePath}/afterBowtieMapping.sh -g ${genomeName} -m ${MERGE_FP} -c ${CONTIG_FP} -d ${DEPTH_FP} -w ${WINDOW} -u ${UNPAIR} -U ${UNFILTER} -p ${PipePath} -M ${MERGE_P} -C ${CONTIG_P} -D ${DEPTH_P} -f ${PLOIDYFILTER} -e ${SINGLE_END} -P ${PEAKCALL} -F ${FOOTPRINT} -W ${WINDOWTRACK} -I ${WINDOWINCR} -b ${BOWTIEMEMORY} 
    
    cd ..
    
echo ""

}

runHubbing(){

#-------Hubbing-starts--------------------------------------------------------------

if [ "${HUBBING}" -eq 1 ] ; then
    
    if [ "${PYRAMIDRERUN}" -eq "0" ]; then
        testedFile="../dnase_pipe_2_param.txt"
        doTempFileTesting
        cp ../dnase_pipe_2_param.txt .
    fi

    # Generating parameter file for the hubbing of this data !
    echo -e "${folderList[$i]}\t${fileList[$i]}\t" >> ./dnase_pipe_1_param.txt
   
#---------Hubbing------------------------------------------------------
    
    printThis=" Hubbing .."
    printToLogFile
    
    cp ../versionInfoHTML.txt .

    # Parameter list defaults..
    hubIsRerun=0
    hubOnlyHub=0
    hubIsPloidyFiltered=1
    

    if [ "${ONLY_FP_AND_PEAK}" -eq 1 ] || [ "${ONLY_PEAK}" -eq 1 ] ; then
        hubIsRerun=1
    fi
    
    if [ "${ONLY_HUB}" -eq 1 ] ; then
        hubOnlyHub=1
    fi
    
    # Here, checking if the run for ploidy filter is actually the PREVIOUS run..
    
    checkPreviousFolder=0
    if [ "${ONLY_FP_AND_PEAK}" -eq 1 ] || [ "${ONLY_PEAK}" -eq 1 ] || [ "${ONLY_HUB}" -eq 1 ] ; then
        checkPreviousFolder=1
    fi
        
    if [ "${checkPreviousFolder}" -eq 0 ]; then
    
    # Normal ploidy filter criteria applies (will be checked from the parameters of THIS run..)
    if [ "${PLOIDYFILTER}" -eq 0 ]; then
        hubIsPloidyFiltered=0
    else
        
        weHavePloidyFile=0
        setPloidyPath
        
        # HERE TO BE REMOVED (IF_CLAUSE) WHEN SUPPORTING OTHER PLOIDY TRACKS TOO !
        if [ "${weHavePloidyFile}" -eq 1 ] ; then    
            hubIsPloidyFiltered=1
        else
            hubIsPloidyFiltered=0
        fi   
    fi
    
    else
    # Normal ploidy filter criteria do not apply - the ploidy filter will be determined from PREVIOUS RUN - grepped out from parameter files saved when originally running filtering..

        # This is rerun for VS14 or higher (we have parameter file saved and fetched)
        if [ -s "parameters_OLD.log" ] ; then
            printThis="Fetching info from earlier run, for hubbing :"
            printToLogFile
            
            printThis="Earlier run was ran in :"
            printToLogFile
            printThis=$( cat parameters_OLD.log | grep VERSION )
            printToLogFile
            
            printThis="Earlier run ploidy tracks parameter :"
            printToLogFile
            printThis=$( cat parameters_OLD.log | grep FILTER | grep -v UNFILTER )
            printToLogFile
            
            # The row in the parameter file is either
            # BLACKLISTFILTER 1  (TRUE=1, FALSE=0)
            # or
            # BLACKLISTFILTER 0  (TRUE=1, FALSE=0)
            hubIsPloidyFiltered=$( cat parameters_OLD.log | grep FILTER | grep -v UNFILTER | sed 's/.*FILTER //' | sed 's/\s(.*$//' )
            
        # This is rerun for VS13 or lower (we don't have parameter file saved and fetched)
        # Now we have only ploidy filtering for mm9, and no flag to "no ploidy filter" - so if genome mm9, we have ploidy filtered, otherwise not.
        else
            if [ "${genomeName}" = "mm9" ]; then
                hubIsPloidyFiltered=1
            else
                hubIsPloidyFiltered=0
            fi
        
        fi
    

    fi
    
    # As the custom genomes need to still be assigned to a "real" genome..
    TEMPORARYgenome="${genomeName}"
    if [ "${TEMPORARYgenome}" == "mm10balb" ]; then
        TEMPORARYgenome="mm10"
    fi

    printThis="${PipePath}/dataHubGenerator.sh -e "${QSUBERRFILE}" -o "${QSUBOUTFILE}" -g "${TEMPORARYgenome}" -n "${MagicNumber}" -w ${WINDOW} -p "${PipePath}" --rerun "${hubIsRerun}" --onlyhub "${hubOnlyHub}" --ploidyTrack "${hubIsPloidyFiltered}" --singleEnd "${SINGLE_END}" --symbolic "${SYMBOLIC}" -W ${WINDOWTRACK} --orangeblue ${ORANGEBLUE} --redgreen ${REDGREEN} "
    printToLogFile

    ${PipePath}/dataHubGenerator.sh -e "${QSUBERRFILE} " -o "${QSUBOUTFILE}" -g "${TEMPORARYgenome}" -n "${MagicNumber}" -w ${WINDOW} -p "${PipePath}" --rerun "${hubIsRerun}" --onlyhub "${hubOnlyHub}" --ploidyTrack "${hubIsPloidyFiltered}" --singleEnd "${SINGLE_END}" --symbolic ${SYMBOLIC} -W ${WINDOWTRACK} --orangeblue ${ORANGEBLUE} --redgreen ${REDGREEN} 
    rm -f versionInfoHTML.txt
    
    rm -f dnase_pipe_1_param.txt
    rm -f dnase_pipe_2_param.txt
  

else
    echo "No hubbing done - as hubbing was not requested ! "
fi
}




runPipe(){
    
# ----------------------------------------------
    
    #-------------------------------------------------
    # This is where BAM and FASTQ pipes merge
    #-------------------------------------------------
    
    # Initial  QC (assumes excisting READ1.fastq, READ2.fastq)
    printThis="Running FastQC for initial set of reads.."
    printToLogFile
    printThis="${PipePath}/QC_and_Trimming.sh --fastqc --single ${SINGLE_END}"
    printToLogFile
    ${PipePath}/QC_and_Trimming.sh --fastqc --single ${SINGLE_END}
    
    # Changing names of fastqc folders to be "ORIGINAL"
    mkdir READ1_fastqc_ORIGINAL
    mv -f READ1_fastqc.html READ1_fastqc_ORIGINAL/fastqc_report.html
    mv -f "READ1_fastqc.zip" "READ1_fastqc_ORIGINAL.zip"
    
    if [ "${SINGLE_END}" -eq 0 ] ; then
    mkdir READ2_fastqc_ORIGINAL
    mv -f READ2_fastqc.html READ2_fastqc_ORIGINAL/fastqc_report.html
    mv -f "READ2_fastqc.zip" "READ2_fastqc_ORIGINAL.zip"
    fi
    
    
    #Check BOWTIE quality scores..
    
    LineCount=$(($( grep -c "" READ1.fastq )/4))
    
    if [ "${LineCount}" -gt 100000 ] ; then
        printThis="bowtieQuals=$( perl ${PerlHelpersPath}/fastq_scores_bowtie${BOWTIE}.pl -i READ1.fastq -r 90000 )"
        printToLogFile
        bowtieQuals=$( perl ${PerlHelpersPath}/fastq_scores_bowtie${BOWTIE}.pl -i READ1.fastq -r 90000 )
    else
        rounds=$((${LineCount}-10))
        printThis="bowtieQuals=$( perl ${PerlHelpersPath}/fastq_scores_bowtie${BOWTIE}.pl -i READ1.fastq -r ${rounds} )"
        printToLogFile
        bowtieQuals=$( perl ${PerlHelpersPath}/fastq_scores_bowtie${BOWTIE}.pl -i READ1.fastq -r ${rounds} )
    fi
    
    echo "Bowtie will be ran in quality score scheme : ${bowtieQuals}"

    # The location of "zero" for the filtering/trimming programs cutadapt, trim_galore, flash    
    intQuals=""
    if [ "${bowtieQuals}" == "--phred33-quals" ] || [ "${bowtieQuals}" == "--phred33" ]; then
        intQuals="33"
    else
        # Both solexa and illumina phred64 have their "zero point" in 64
        intQuals="64"
    fi


    #Run BOWTIE to all builds requested, before deleting fastqs :
    for j in $( seq 0 $((${#GENOMEARRAY[@]} - 1)) ); do
    #for j in "${GENOMEARRAY[@]}"; do
        genomeName=${GENOMEARRAY[$j]}
        
        setBOWTIEgenomeSizes
        runBowtie
        echo
        
    done

    #Delete FASTQs (after checking that bowtie succeeded)
    deleteFastq
    echo
    
    #Run the rest of the pipe
    for j in $( seq 0 $((${#GENOMEARRAY[@]} - 1)) ); do
    #for j in "${GENOMEARRAY[@]}"; do
        genomeName=${GENOMEARRAY[$j]}
        
        runBowtieSubscript
        
        runHubbing
        
        cd ${genomeName}
        cleanUpFolder
        listFolder
        echo
        cat hub_address.txt
        cp ../../parameters.log ForPeakCallAndFootPrintRerun/.
        cd ..
        
        echo
    done
    
}


# SCRIPT OUTLINE

# Script makes bam --> fastq transition (if starting with bam files),
# or fetches fastq files (if starting with fastq files).
# 
# After that runs bowtie.
# Unmapping reads continue to trimming, and short-read combining via flash.
# Generates fastqc folders before and after trimming/flashing, an for unmapping reads.
# Unmapping reads are then mapped as single end reads (if paired end sample).
# 
# Generates filtered.bam (no ploidy, no duplicates, only proper pairs), returns it to original fastq input order.
# For both raw.bam and filtered.bam regenerates fragments (only proper pairs, if paired sample), generates bedgraph pileup, and its bigwig transform.
#
# For filtered.bam generates peak calls, and DNase footprints (leftmost 1base LEFT reads, rightmost 1base RIGHT reads) - in bedgraph and bigwig formats.
#
# Generates data UCSC data hub including i) the bigwigs, ii) the bowtie analysis statistics, and iii) the fastqc folders.
#

# PARAMETER FILE FORMAT IN DETAIL

# The script demands parameter files.
# The parameter files should be <TAB> delimited, and end with a SINGLE EMPTY LINE, and be located in the running folder
# Other whitespace is permitted in addition to <TAB> , and multiple consequtive whitespace are collapsed into one.

# PIPE_bamPaths.txt lists the sample names, and sample paths (bam files).
#
# MyCellLine_DNaseI 	/hts/data3/telenius/sample_first.bam
# OtherCellLine_H3K4me1 	/hts/data6/otherUser/otherFolder/sample_second_file.bam
#

# PIPE_fastqPaths.txt lists the sample names, and sample paths (fastqc files).
#
# Format 1 ("old format" - one lane per sample)
# MyCellLine_DNaseI 	/hts/data3/telenius/sample_first_R1.fastq   /hts/data3/telenius/sample_first_R2.fastq
# OtherCellLine_H3K4me1 	/hts/data6/otherUser/otherFolder/sample_second_file_READ1.fastq    /hts/data6/otherUser/otherFolder/sample_second_file_READ2.fastq
#
# Format 1 ("old format" - many lanes per sample)
# MyCellLine_DNaseI 	/hts/data3/telenius/sample_first_R1_L1.fastq,/hts/data3/telenius/sample_first_R1_L2.fastq   /hts/data3/telenius/sample_first_R2_L1.fastq,/hts/data3/telenius/sample_first_R2_L2.fastq
#
#
# Format 2 ("new format" - one lane per sample)
# MyCellLine_DNaseI 	sample_first_R1.fastq   sample_first_R2.fastq      /hts/data3/telenius
# OtherCellLine_H3K4me1 	sample_second_file_READ1.fastq  sample_second_file_READ2.fastq   /hts/data6/otherUser/otherFolder
#
# Format 2 ("new format" - many lanes per sample)
# MyCellLine_DNaseI 	sample_first_R1_L1.fastq,sample_first_R1_L2.fastq   sample_first_R2_L1.fastq,sample_first_R2_L2.fastq   /hts/data3/telenius
#

# PIPE_hubbingSymbolic.txt defines the public folders for data visualisation, and hub generation (oneliner, ending with newline) :
#
# MyHubName /public/telenius/MyNewDataHub  /permanent/disk/location/for/BigWigs
#
# This hub will be called "MyHubName" in UCSC , and it will be located in servertype://server.address.ex.amp.le/public/telenius/MyNewDataHub/hub.txt
# See the conf/config.sh to set up the editing of the second column here : to match your own server name and folder structure in your public disk.
# The location you give here is the DISK PATH to the location (not the server address - if they should differ in the file path)
#
# If the folder in column 2 already exists, it is assumed this is a RERUN (peak call, footprint), and results are APPLIED to existing data hub.
# The only exception are the description.html and tracks.txt files, which are overwritten BUT old file is safety-copied to the folder before overwrite.
# If the folder in column 2 does NOT exist, a new data hub is generated from scratch.
#
# The bigwigs generated for the hub will be stored in folder /permanent/disk/location/for/BigWigs/MyHubName
# This location should be visible to the public disk space (column 2) - as symbolic links are created to make the stored files "appear to be" in the public server space.
# This disk location should be such, that it will NOT be moved around (hubs will break if these files are moved elsewhere - as the mentioned symbolic links will break).
# It is recommended to ALWAYS use the same disk location as the 3rd column of the file (to easier to remember not to move this folder around).
# If folder /permanent/disk/location/for/BigWigs/MyHubName already exists, the files within it are overwritten (if they have same names as existing files)
#

# PIPE_hubbing.txt defines the public folders for data visualisation, and hub generation (oneliner, ending with newline) :
#
# MyHubName /public/telenius/MyNewDataHub
#
# This data hub is exactly the same as the above, but does NOT generate symbolic links to the bigwig files.
# Instead - stores the bigwigs in the public disk storage location (given in column 2).
# This is not recommended, as the bigwig files are reasonably large, and thus small public area disks may get full quickly.
#
#

######################################################################
# The codes of the pipeline 
######################################################################
#
# NGseqBasic
#
# |-- NGseqBasic.sh
# |
# `-- bin
#     |-- mainScripts
    
#     |   |-- QC_and_Trimming.sh
#     |   |-- afterBowtieMapping.sh 
#     |   `-- dataHubGenerator.sh
#     |   
#     |-- bashHelpers
#     |   |-- inputFastqs.sh
#     |   |-- mappingSubroutines.sh
#     |   |-- trimmingSubroutines.sh
#     |   |-- bowtieStatistics.sh (used by mappingSubroutines.sh )
#     |   |-- cleanUpAndList.sh (used by NGseqBasic.sh)
#     |   |-- fileTesters.sh (used by most scripts)
#     |   |-- logFilePrinter.sh (used by most scripts)
#     |   |-- parameterSetters.sh (used by NGseqBasic.sh)
#     |   `-- usageAndVersion.sh (used by NGseqBasic.sh)
#     |   
#     `-- perlHelpers
#         |-- data2gff.pl (used by afterBowtieMapping.sh )
#         |-- fastq_scores_bowtie1.pl (used by trimmingSubroutines.sh)
#         |-- fastq_scores_bowtie2.pl (used by trimmingSubroutines.sh)
#         |-- reverse_seq.pl (used by trimmingSubroutines.sh)
#         |-- sam2fastq.pl (used by NGseqBasic.sh)
#         |-- trim1base3prime.pl (used by trimmingSubroutines.sh)
#         `-- windowingScript.pl (used by afterBowtieMapping.sh)
#
#
######################################################################

QSUBOUTFILE="qsub.out"
QSUBERRFILE="qsub.err"

LANES=1
BOWTIE=1
SKIP_BOWTIE=0
FOOTPRINT=0
PEAKCALL=0
ONLY_HUB=0
ONLY_FP_AND_PEAK=0
ONLY_PEAK=0
# Possible values for these are $(( 00 ))  = neither, $(( 01 )) = only mapped, $(( 10 )) = all reads, $(( 11 )) = both
UNFILTER=0
UNTRIMMED=0
# Possible values for this are 00 = neither, 01 = filtered, 10 = unfiltered, 11 = both
UNPAIR=0
SAVEUNMAPPED=0
TRIM=1
FLASH=1
GENOMES=""
NORMDEPTH=0
MAXINS=350
MERGE_FP=100
CONTIG_FP=10
DEPTH_FP=10
MERGE_P=0
CONTIG_P=20
DEPTH_P=40
WINDOW=1
DOWINDOWING=1
WINDOWTRACK=0
WINDOWINCR=0
CAPITAL_M=0
LOWERCASE_M=0
BOWTIEMEMORY=256
BOWTIE_PROCESSORS=3
otherBowtie1Parameters=""
otherBowtie2Parameters=""
LOWERCASE_V=-1
bowtie1MismatchBehavior=""
bowtie2MismatchBehavior=""
PLOIDYFILTER=1
SINGLE_END=0
SAVE_BDG_AND_WIG=0
GZIP=0
TRIM5=0
NEXTERA=0
ORANGEBLUE=0
REDGREEN=0

PYRAMIDRERUN=0

CUSTOMAD=-1
ADA31="no"
ADA32="no"
ADA51="no"
ADA52="no"

genomeName=""
genomeBuild=""
bowtieGenomeBuild=""

GENOMEARRAY[0]="UNDEFINED"


# This binary 1/0 tells if we originated to RunPipe from fastq or bam parameter file (and determines if the raw bam file is saved or not
# - if we start from bams, we already HAVE it and thus do not save)
fastqBasedRun=0

dirForQuotaAsking=$( pwd )

#------------------------------------------

# Loading subroutines in ..

PipeTopPath="$( which $0 | sed 's/\/NGseqBasic.sh$//' )"

BashHelpersPath="${PipeTopPath}/bin/bashHelpers"

# CLEANING UP AND LISTING FOLDER CONTENTS
. ${BashHelpersPath}/cleanUpAndList.sh
# TESTING FILE EXISTENCE
. ${BashHelpersPath}/fileTesters.sh
# PRINTING TO LOG AND ERROR FILES
. ${BashHelpersPath}/logFilePrinter.sh
# SETTING THE RUN PARAMETERS
. ${BashHelpersPath}/parameterSetters.sh
# SETTING THE GENOME BUILD PARAMETERS
. ${BashHelpersPath}/genomeSetters.sh
# SETTING THE BLACKLIST GENOME LIST PARAMETERS
. ${BashHelpersPath}/blacklistSetters.sh
# PRINTING HELP AND VERSION MESSAGES
. ${BashHelpersPath}/usageAndVersion.sh

# LOADING FASTQS AND COMBINING LANES
. ${BashHelpersPath}/inputFastqs.sh
# MAPPING IN BOWTIE
. ${BashHelpersPath}/mappingSubroutines.sh
# TRIMMING SUBROUTINES
. ${BashHelpersPath}/trimmingSubroutines.sh

#------------------------------------------

# From where to call the main scripts operating from the mainScripts folder..

PipePath="${PipeTopPath}/bin/mainScripts"

#------------------------------------------

# From where to call the helper perl scripts..

PerlHelpersPath="${PipeTopPath}/bin/perlHelpers"

#------------------------------------------

# From where to call the CONFIGURATION script..

confFolder="${PipeTopPath}/conf"

#------------------------------------------

# Calling in the CONFIGURATION script and its default setup :

supportedGenomes=()
BOWTIE1=()
BOWTIE2=()
UCSC=()
genomesWhichHaveBlacklist=()
BLACKLIST=()

. ${confFolder}/loadNeededTools.sh
. ${confFolder}/genomeBuildSetup.sh
. ${confFolder}/serverAddressAndPublicDiskSetup.sh

setGenomeLocations
setPublicLocations

#------------------------------------------


OPTS=`getopt -o h,m:,M:,p:,v: --long help,test,veryOldChIP,saveFastq,orangeBlue,redGreen,gz,singleEnd,footPrint,peakCall,nextera,bowtie1,bowtie2,seedmms:,seedlen:,maqerr:,noBowtie,blacklistFilter,noBlacklistFilter,onlyHub,onlyPeakCall,onlyFPandPC,saveUnmapped,saveUnpaired,saveBDG,saveUnpairedFiltered,saveUnfiltered,saveUnfilteredMapped,saveUntrimmed,saveUntrimmedMapped,trim,noTrim,flash,noFlash,normDepth,noWindow,trim5,outfile:,errfile:,chunkmb:,lanes:,genomes:,maxins:,mergeFP:,contigFP:,depthFP:,mergeP:,contigP:,depthP:,windowSize:,windowIncr:,windowFP:,ada3read1:,ada3read2:,ada5read1:,ada5read2:,pyramidRerun -- "$@"`
if [ $? != 0 ]
then
    usage ;
fi

eval set -- "$OPTS"

while true ; do
    case "$1" in
        -h) usage ; shift;;
        -m) LOWERCASE_M=$2 ; shift 2;;
        -M) CAPITAL_M=$2 ; shift 2;;
        -p) BOWTIE_PROCESSORS=$2 ; shift 2;;
        -v) LOWERCASE_V=$2; shift 2;;
        --help) usage ; shift;;
        --test) echo "TEST RUN - PYRAMID DEVELOPMENT !"; shift;;
        --veryOldChIP) echo "Pyramid run for ANCIENT data - no knowledge of antibody supplier or batch id."; shift;;
        --saveFastq) echo "Pyramid run for GEO data, where fastq files were also saved."; shift;;
        --singleEnd) SINGLE_END=1 ; shift;;
        --orangeBlue) ORANGEBLUE=1 ; shift;;
        --redGreen) REDGREEN=1 ; shift;;
        --chunkmb) BOWTIEMEMORY=$2 ; shift 2;;
        --lanes) LANES=$2 ; shift 2;;
        --nextera) NEXTERA=1 ; shift;;
        --bowtie1) BOWTIE=1 ; shift;;
        --bowtie2) BOWTIE=2 ; shift;;
        --noBowtie) SKIP_BOWTIE=1 ; shift;;
        --noBlacklistFilter) PLOIDYFILTER=0 ; shift;;
        --blacklistFilter) PLOIDYFILTER=1 ; shift;;
        --footPrint) FOOTPRINT=1; shift;;
        --peakCall) PEAKCALL=1; shift;;
        --pyramidRerun) PYRAMIDRERUN=1; shift;;
        --onlyHub) ONLY_HUB=1 ; shift;;
        --onlyPeakCall) ONLY_PEAK=1 ; shift;;
        --onlyFPandPC) ONLY_FP_AND_PEAK=1 ; shift;;
        --saveBDG) SAVE_BDG_AND_WIG=1 ; shift;;
        --saveUnmapped) SAVEUNMAPPED=1 ; shift;;
        --saveUnpaired) UNPAIR=$(( ${UNPAIR}+10 )) ; shift;;
        --saveUnpairedFiltered) UNPAIR=$(( ${UNPAIR}+1 )) ; shift;;
        --saveUnfiltered) UNFILTER=$(( ${UNFILTER}+10 )) ; shift;;
        --saveUnfilteredMapped) UNFILTER=$(( ${UNFILTER}+1 )) ; shift;;
        --saveUntrimmed) UNTRIMMED=$(( ${UNTRIMMED}+10 )) ; shift;;
        --saveUntrimmedMapped) UNTRIMMED=$(( ${UNTRIMMED}+1 )) ; shift;;
        --trim) TRIM=1 ; shift;;
        --noTrim) TRIM=0 ; shift;;
        --trim5) TRIM5=1 ; shift;;
        --flash) FLASH=1 ; shift;;
        --noFlash) FLASH=0 ; shift;;
        --noWindow) DOWINDOWING=0 ; shift;;
        --outfile) QSUBOUTFILE=$2 ; shift 2;;
        --errfile) QSUBERRFILE=$2 ; shift 2;;
        --windowIncr) WINDOWINCR=$2 ; shift 2;;
        --windowSize) WINDOWTRACK=$2 ; shift 2;;
        --normDepth) NORMDEPTH=1 ; shift;;
        --genomes) GENOMEARRAY=($(echo $2 | sed 's/,/\n/g')) ; shift 2;;
        --maxins) MAXINS=$2 ; shift 2;;
        --mergeFP) MERGE_FP=$2 ; shift 2;;
        --contigFP) CONTIG_FP=$2 ; shift 2;;
        --depthFP) DEPTH_FP=$2 ; shift 2;;
        --mergeP) MERGE_P=$2 ; shift 2;;
        --contigP) CONTIG_P=$2 ; shift 2;;
        --depthP) DEPTH_P=$2 ; shift 2;;
        --windowFP) WINDOW=$2 ; shift 2;;
        --incr) INCREMENT=$2 ; shift 2;;
        --gz) GZIP=1 ; shift;;
        --ada3read1) ADA31=$2 ; shift 2;;
        --ada3read2) ADA32=$2 ; shift 2;;
        --ada5read1) ADA51=$2 ; shift 2;;
        --ada5read2) ADA52=$2 ; shift 2;;
        --seedmms) bowtie1MismatchBehavior="${bowtie1MismatchBehavior} --seedmms $2 " ; ${bowtie2MismatchBehavior}="${bowtie2MismatchBehavior} -N $2 "  ; shift 2;;
        --seedlen) bowtie1MismatchBehavior="${bowtie1MismatchBehavior} --seedlen $2 " ; ${bowtie2MismatchBehavior}="${bowtie2MismatchBehavior} -L $2 " ; shift 2;;
        --maqerr) bowtie1MismatchBehavior="${bowtie1MismatchBehavior} --maqerr $2 " ; shift 2;;
        --) shift; break;;
    esac
done

#------------------------------------------

echo "NGseqBasic.sh - by Jelena Telenius, 03/11/2016"
echo
timepoint=$( date )
echo "run started : ${timepoint}"
echo
echo "Script located at"
which $0
echo

echo "RUNNING IN MACHINE : "
hostname --long

echo "run called with parameters :"
echo "NGseqBasic.sh" $@
echo

#------------------------------------------

echo
echo "PipeTopPath ${PipeTopPath}"
echo "PipePath ${PipePath}"
echo "confFolder ${confFolder}"
echo "BashHelpersPath ${BashHelpersPath}"
echo "PerlHelpersPath ${PerlHelpersPath}"
echo

# ------------------------------------------

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

# #####################################################
# 
# NOTE !! about for-loops ( within NGseqBasic script and the scripts sharing its name space ).
#
#
# MASTER LOOPS :
#
# Outerloop = i ( folderList fileList )
# Outerloop = j ( genomelist , if needed = for bam and fastq files, which can map to multiple genomes)
#
# Do not use $i or $j for ANYTHING else - these are the "main loops" and not to be disturbed !
#
# PARAMETER LISTING (and other secondary importance "temporary" loops) :
#
# Lanes = l
# Genomes = g
# Files = f
# Other = k
# 
#
# #####################################################


# ----------------------------------------------

# Modifying and adjusting parameter values, based on run flags

setParameters

# ----------------------------------------------

# Loading the environment - either with module system or setting them into path.
# This subroutine comes from conf/config.sh file

printThis="LOADING RUNNING ENVIRONMENT"
printToLogFile

setPathsForPipe


#---------------------------------------------------------

# If no genomes given (the only obligatory parameter) - if this is not a run for already-mapped data

if [ ${GENOMEARRAY[0]} == "UNDEFINED" ] && [ ! -r "./PIPE_mappedBamPaths.txt" ] && [ ! -r "./PIPE_previousRunPaths.txt" ] ; then

    echo  >&2
    echo "No genomes given with parameter --genomes in the run.sh file - pipeline aborted"  >&2
    usage
    echo  >&2
    exit 1

fi

# ----------------------------------------------

MagicNumber=$$



#--------Generating-the-parameter-files-for-the-subscripts------------------------------------------------------

if [ ! -r "./PIPE_bamPaths.txt" ] && [ ! -r "./PIPE_fastqPaths.txt" ] && [ "${SKIP_BOWTIE}" -eq 0 ] && [ "${ONLY_HUB}" -eq 0 ] && [ "${ONLY_PEAK}" -eq 0 ] && [ "${ONLY_FP_AND_PEAK}" -eq 0 ]; then
    echo  >&2
    echo "No parameter file PIPE_bamPaths.txt or PIPE_fastqPaths.txt provided - pipeline aborted"  >&2
        
    echo  >&2
    exit 1 
elif [ ! -r "./PIPE_mappedBamPaths.txt" ] && [ "${SKIP_BOWTIE}" -eq 1 ] ; then
    echo  >&2
    echo "No parameter file PIPE_mappedBamPaths.txt (for existing bowtie-aligned bam paths) provided - pipeline aborted"  >&2
    parameterUsage
    echo  >&2
    exit 1
elif [ ! -r "./PIPE_previousRunPaths.txt" ] && [ "${ONLY_PEAK}" -eq 1 ] ; then
    echo  >&2
    echo "No parameter file PIPE_previousRunPaths.txt (for existing bowtie-aligned bam paths) provided - pipeline aborted"  >&2
    parameterUsage
    echo  >&2
    exit 1
elif [ ! -r "./PIPE_previousRunPaths.txt" ] && [ "${ONLY_FP_AND_PEAK}" -eq 1 ] ; then
    echo  >&2
    echo "No parameter file PIPE_previousRunPaths.txt (for existing bowtie-aligned bam paths) provided - pipeline aborted"  >&2
    parameterUsage
    echo  >&2
    exit 1
elif [ ! -r "./PIPE_previousRunPaths.txt" ] && [ "${ONLY_HUB}" -eq 1 ] ; then
    echo  >&2
    echo "No parameter file PIPE_previousRunPaths.txt (for existing bowtie-aligned bam paths) provided - pipeline aborted"  >&2
    parameterUsage
    echo  >&2
    exit 1
elif [ ! -r "./PIPE_hubbing.txt" ] && [ ! -r "./PIPE_hubbingSymbolic.txt" ] && [ "${ONLY_HUB}" -eq 1 ] ; then
    echo  >&2
    echo "No parameter file PIPE_hubbing.txt or PIPE_hubbingSymbolic.txt provided - pipeline aborted"  >&2
    parameterUsage
    echo  >&2
    exit 1
fi    

#---------------------------------------------------------
# Here parsing the parameter files - if they are not purely tab-limited, but partially space-limited, or multiple-tab limited, this fixes it.

echo
echo "PARAMETER FILES GIVEN IN RUN FOLDER :"
echo

for file in ./PIPE*.txt
    do
        echo ${file}
        sed -i 's/\s\s*/\t/g' ${file}
    done
    
#------------------------------------------------------
# Here reading hubbing parameters in..

HUBBING=0

if [ -r "./PIPE_hubbing.txt" ] ; then
    cp -f ./PIPE_hubbing.txt ./dnase_pipe_2_param.txt
    HUBBING=1
fi

SYMBOLIC=0

if [ -r "./PIPE_hubbingSymbolic.txt" ] ; then
    cp -f ./PIPE_hubbingSymbolic.txt ./dnase_pipe_2_param.txt
    HUBBING=1
    SYMBOLIC=1
fi

echo
if [ "${HUBBING}" -eq 1 ]; then
    echo "Data hub generation is requested."
    
    if [ "${SYMBOLIC}" -eq 1 ]; then
    echo "Will use symbolic links to store data in public area."
    fi
else
    echo "Data hub generation is NOT requested."
    echo "No UCSC data hub will be generated !"    
fi
echo

#------------------------------------------------------


echo
echo "Run with parameters :"
echo ""

echo "Output log file ${QSUBOUTFILE}" > parameters.log
echo "Output error log file ${QSUBERRFILE}" >> parameters.log
echo "fastqBasedRun ${fastqBasedRun}  (TRUE=1, FALSE=0)" >> parameters.log
echo "SINGLE_END ${SINGLE_END}  (TRUE=1, FALSE=0)" >> parameters.log
echo "LANES ${LANES}" >> parameters.log
echo "GZIP ${GZIP} (TRUE=1, FALSE=0) - if fastq input files are gzipped "
echo ""  >> parameters.log
echo "SKIP_BOWTIE ${SKIP_BOWTIE}  (TRUE=1, FALSE=0)" >> parameters.log
echo "BOWTIE ${BOWTIE} (Bowtie1=1, Bowtie2=2)" >> parameters.log
echo "BOWTIE_PROCESSORS ${BOWTIE_PROCESSORS}" >> parameters.log
echo "BOWTIEMEMORY ${BOWTIEMEMORY}" >> parameters.log
echo "MAXINS ${MAXINS}" >> parameters.log
echo "CAPITAL_M ${CAPITAL_M}"  >> parameters.log
echo "LOWERCASE_M ${LOWERCASE_M}" >> parameters.log
if [ "${BOWTIE}" -eq 1 ]; then
echo "otherBowtie1Parameters ${otherBowtie1Parameters}" >> parameters.log
else
echo "otherBowtie2Parameters ${otherBowtie2Parameters}" >> parameters.log
fi
echo ""  >> parameters.log
echo "TRIM ${TRIM}  (TRUE=1, FALSE=0)" >> parameters.log
echo "TRIM5 ${TRIM5}  (TRUE=1, FALSE=0)" >> parameters.log
echo "NEXTERA ${NEXTERA} (TRUE=1, FALSE=0)" >> parameters.log
echo "CUSTOMAD ${CUSTOMAD}" >> parameters.log

if [ "${CUSTOMAD}" -ne -1 ]; then

echo "ADA31 ${ADA31}" >> parameters.log
echo "ADA32 ${ADA32}" >> parameters.log
echo "ADA51 ${ADA51}" >> parameters.log
echo "ADA52 ${ADA52}" >> parameters.log
   
fi

echo "FLASH ${FLASH}  (TRUE=1, FALSE=0)" >> parameters.log
echo "BLACKLISTFILTER ${PLOIDYFILTER}  (TRUE=1, FALSE=0)" >> parameters.log
echo ""  >> parameters.log
echo "SAVE_BDG_AND_WIG ${SAVE_BDG_AND_WIG}  (TRUE=1, FALSE=0)" >> parameters.log
echo "UNMAPPED ${UNMAPPED}  (TRUE=1, FALSE=0) - save FASTQ of unmapped (SE) / not-mapped-in-proper-pair (PE)" >> parameters.log
echo "UNPAIR ${UNPAIR}  (FALSE=0, filteredReads=1 unfilteredReads=10, both=11) " >> parameters.log
echo "UNFILTER ${UNFILTER}  (FALSE=0, mappedReads=1 allReads=10, both=11) - save UNFILTERED bam" >> parameters.log
echo "UNTRIMMED ${UNTRIMMED}  (FALSE=0, mappedReads=1 allReads=10, both=11) - save UNTRIMMED bam" >> parameters.log
echo ""  >> parameters.log
echo "DOWINDOWING ${DOWINDOWING}  (TRUE=1, FALSE=0) - windowing for the red and green graphs" >> parameters.log
echo "WINDOWTRACK ${WINDOWTRACK}" >> parameters.log
echo "WINDOWINCR ${WINDOWINCR}" >> parameters.log
echo ""  >> parameters.log
echo "FOOTPRINT ${FOOTPRINT}  (TRUE=1, FALSE=0)" >> parameters.log
echo "MERGE_FP ${MERGE_FP}" >> parameters.log
echo "CONTIG_FP ${CONTIG_FP}" >> parameters.log
echo "DEPTH_FP ${DEPTH_FP}" >> parameters.log
echo "WINDOW ${WINDOW}" >> parameters.log
echo ""  >> parameters.log
echo "PEAKCALL ${PEAKCALL}  (TRUE=1, FALSE=0)" >> parameters.log
echo "MERGE_P ${MERGE_P}" >> parameters.log
echo "CONTIG_P ${CONTIG_P}" >> parameters.log
echo "DEPTH_P ${DEPTH_P}" >> parameters.log
echo ""  >> parameters.log
echo "PYRAMIDRERUN ${PYRAMIDRERUN}  (TRUE=1, FALSE=0)" >> parameters.log
echo "ONLY_FP_AND_PEAK ${ONLY_FP_AND_PEAK} (TRUE=1, FALSE=0)" >> parameters.log
echo "ONLY_PEAK ${ONLY_PEAK}  (TRUE=1, FALSE=0)" >> parameters.log
echo "ONLY_HUB ${ONLY_HUB}  (TRUE=1, FALSE=0)" >> parameters.log
echo "HUBBING ${HUBBING}  (TRUE=1, FALSE=0)" >> parameters.log
echo "SYMBOLIC ${SYMBOLIC}" >> parameters.log
echo "ORANGEBLUE ${ORANGEBLUE}" >> parameters.log
echo "REDGREEN ${REDGREEN} (the old non-colorblind-friendly default red-green colors)" >> parameters.log
echo ""  >> parameters.log

# The public paths are just listed here, not used.
# They will be separately fetched again, in the dataHubGenerator.sh
echo "SERVERTYPE ${SERVERTYPE}" >> parameters.log
echo "SERVERADDRESS ${SERVERADDRESS}" >> parameters.log
echo "ADDtoPUBLICFILEPATH ${ADDtoPUBLICFILEPATH}" >> parameters.log
echo "REMOVEfromPUBLICFILEPATH ${REMOVEfromPUBLICFILEPATH}" >> parameters.log
echo "tobeREPLACEDinPUBLICFILEPATH ${tobeREPLACEDinPUBLICFILEPATH}" >> parameters.log
echo "REPLACEwithThisInPUBLICFILEPATH ${REPLACEwithThisInPUBLICFILEPATH}" >> parameters.log

cat parameters.log
echo

version
echo -e ${versionInfo} >> parameters.log


if [ "${TRIM}" -eq 0 ] && [ "${FLASH}" -eq 1 ]; then
    echo "NOTE !! INTERPRET YOUR RESULTS WITH CAUTION !! :"
    echo "Combination --flash --noTrim is NOT RECOMMENDED (can result in combining reads on the sites of ADAPTERS instead of the reads themselves). "
    echo
fi


#--------This-run-is-a-RERUN-either-to-optimise-the-parameters-or-simply-create-DATA-HUB-----------------

if [ "${ONLY_FP_AND_PEAK}" -eq 1 ] || [ "${ONLY_PEAK}" -eq 1 ] || [ "${ONLY_HUB}" -eq 1 ] ; then

folderList=($( cut -f 1 ./PIPE_previousRunPaths.txt ))
fileList=($( cut -f 2 ./PIPE_previousRunPaths.txt ))
GENOMEARRAY=($( cut -f 3 ./PIPE_previousRunPaths.txt ))

printRunStartArrays
    
        #Run the rest of the pipe
        for i in $( seq 0 $((${#GENOMEARRAY[@]} - 1)) ); do
        #for i in "${GENOMEARRAY[@]}"; do
            genomeName=${GENOMEARRAY[i]}
            
            printThis="Starting run for sample : ${folderList[$i]}, mapped to ${genomeName}"
            printToLogFile
            
            TRUEgenomeName="${genomeName}"
            # Here take into account that for pyramid runs we only have one genome, and we don't have a folder for it..
            if [ "${PYRAMIDRERUN}" -eq "1" ]; then
                genomeName=""   
            fi
            
            mkdir ${folderList[$i]}
            mkdir ${folderList[$i]}/${TRUEgenomeName}
            
            cd ${folderList[$i]}/${TRUEgenomeName}
            pwd
            
            # The existence of "ForPeakCallAndFootPrintRerun" folder is the TRIGGER of this part of the script in the bowtiesubscript !
            mkdir Rerun
            
            #-RERUN-FOR-BOTH-FOOTPRINT-AND-PEAK-CALL-------
            if [ "${ONLY_FP_AND_PEAK}" -eq 1 ] ; then
                
                # Backwards compatible - if BED files exist, using them. Otherwise making bed files from bigbeds.
                
                if  [ -s ${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/RIGHT_sortedTempFile.bed ]; then
                    cp -f ${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/RIGHT_sortedTempFile.bed Rerun/RIGHT_sortedTempFile.bed
                else
                    cd Rerun
                    bigBedToBed ${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/RIGHT_sortedTempFile.bb RIGHT_sortedTempFile.bed
                    cd ..
                fi
                
                if  [ -s ${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/LEFT_sortedTempFile.bed ]; then
                    cp -f ${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/LEFT_sortedTempFile.bed Rerun/LEFT_sortedTempFile.bed
                else
                    cd Rerun
                    bigBedToBed ${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/LEFT_sortedTempFile.bb LEFT_sortedTempFile.bed
                    cd ..
                fi
                
                if  [ -s ${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/READ_sortedTempFile.bed ]; then
                    cp -f ${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/READ_sortedTempFile.bed Rerun/READ_sortedTempFile.bed
                else
                    cd Rerun
                    bigBedToBed ${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/READ_sortedTempFile.bb READ_sortedTempFile.bed
                    cd ..
                fi
                
                cp -f ${fileList[$i]}/${genomeName}/BigWigs/filtered_pileup_scaled.bw .
                
                echo "Fetched files :"
                echo "${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/RIGHT_sortedTempFile.bed "
                echo "${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/LEFT_sortedTempFile.bed "
                echo "${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/READ_sortedTempFile.bed "
                
        # -g) genomeName=$2 ; shift 2;;
        # -m) MERGE_FP=$2 ; shift 2;;
        # -c) CONTIG_FP=$2 ; shift 2;;
        # -d) DEPTH_FP=$2 ; shift 2;;
        # -w) WIN=$2 ; shift 2;;
        # -u) saveUnpaired=$2 ; shift 2;;
        # -U) saveUnfiltered=$2 ; shift 2;;
        # -p) PipePath=$2 ; shift 2;;
        # -M) MERGE_P=$2 ; shift 2;;
        # -C) CONTIG_P=$2 ; shift 2;;
        # -D) DEPTH_P=$2 ; shift 2;;
        # -f) PLOIDYFILTER=$2 ; shift 2;;
        # -e) SINGLEEND=$2 ; shift 2;;
        # -P) PEAKCALL=$2 ; shift 2;;
        # -F) FOOTPRINT=$2 ; shift 2;;
        # -W) WINTRACK=$2 ; shift 2;;
                
                printThis="${PipePath}/afterBowtieMapping.sh -g ${TRUEgenomeName} -m ${MERGE_FP} -c ${CONTIG_FP} -d ${DEPTH_FP} -w ${WINDOW} -u ${UNPAIR} -U ${UNFILTER} -p ${PipePath} -M ${MERGE_P} -C ${CONTIG_P} -D ${DEPTH_P} -f ${PLOIDYFILTER} -e ${SINGLE_END} -P ${PEAKCALL} -F ${FOOTPRINT} -W ${WINDOWTRACK} -I ${WINDOWINCR}"
                printToLogFile
                ${PipePath}/afterBowtieMapping.sh -g ${TRUEgenomeName} -m ${MERGE_FP} -c ${CONTIG_FP} -d ${DEPTH_FP} -w ${WINDOW} -u ${UNPAIR} -U ${UNFILTER} -p ${PipePath} -M ${MERGE_P} -C ${CONTIG_P} -D ${DEPTH_P} -f ${PLOIDYFILTER} -e ${SINGLE_END} -P ${PEAKCALL} -F ${FOOTPRINT} -W ${WINDOWTRACK} -I ${WINDOWINCR}

            fi
            
            #-RERUN-FOR-PEAK-CALL--------------------------
            if [ "${ONLY_PEAK}" -eq 1 ]  ; then
                
                if  [ -s ${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/READ_sortedTempFile.bed ]; then
                    cp -f ${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/READ_sortedTempFile.bed Rerun/READ_sortedTempFile.bed
                else
                    cd Rerun
                    bigBedToBed ${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/READ_sortedTempFile.bb READ_sortedTempFile.bed
                    cd ..
                fi
                
                echo "Fetched files :"
                echo "${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/READ_sortedTempFile.bed "
                
                 printThis="${PipePath}/afterBowtieMapping.sh -g ${TRUEgenomeName} -m ${MERGE_FP} -c ${CONTIG_FP} -d ${DEPTH_FP} -w ${WINDOW} -u ${UNPAIR} -U ${UNFILTER} -p ${PipePath} -M ${MERGE_P} -C ${CONTIG_P} -D ${DEPTH_P} -f ${PLOIDYFILTER} -e ${SINGLE_END} -P ${PEAKCALL} -F ${FOOTPRINT} -W ${WINDOWTRACK} -I ${WINDOWINCR}"
                printToLogFile
                ${PipePath}/afterBowtieMapping.sh -g ${TRUEgenomeName} -m ${MERGE_FP} -c ${CONTIG_FP} -d ${DEPTH_FP} -w ${WINDOW} -u ${UNPAIR} -U ${UNFILTER} -p ${PipePath} -M ${MERGE_P} -C ${CONTIG_P} -D ${DEPTH_P} -f ${PLOIDYFILTER} -e ${SINGLE_END} -P ${PEAKCALL} -F ${FOOTPRINT} -W ${WINDOWTRACK} -I ${WINDOWINCR}
            fi
            
            rm -rf Rerun
            
            #-GENERATION-OF-DATA-HUB------------------------
            if [ "${ONLY_HUB}" -eq 1 ]  ; then
                
                WINDOW="NA"
                    
                #ln -fs "${fileList[$i]}/${genomeName}/BigWigs/*.bw" . 
                #ln -fs  "${fileList[$i]}/READ*_fastqc_ORIGINAL" ../.
                #ln -fs  "${fileList[$i]}/READ*_fastqc_TRIMMED" ../.
                #ln -fs "${fileList[$i]}/read_trimming.log" ../.
                #ln -fs "${fileList[$i]}/${genomeName}/SamTools_statistics/*_statistics.log" .
                
                for filename in ${fileList[$i]}/${genomeName}/BigWigs/*.bw
                do
                    #echo ${filename}
                    cp -f ${filename} .   
                done
                
                for filename in ${fileList[$i]}/${genomeName}/SamTools_statistics/*_statistics.log
                do
                    #echo ${filename}
                    cp -f ${filename} .   
                done
                
                for filename in ${fileList[$i]}/READ*_fastqc_ORIGINAL
                do
                    #echo ${filename}
                    cp -rf ${filename} ../.   
                done
                
                for filename in ${fileList[$i]}/READ*_fastqc_TRIMMED
                do
                    #echo ${filename}
                    cp -rf ${filename} ../.
                done
                
                cp -f ${fileList[$i]}/read_trimming.log ../.
                
                
                echo "Copied files (in genome folder):"
                ls -l
                echo "Copied files (in upper folder):"
                ls -l ../.
                
                cd ..
                
                # For hubbing - info file copying.. (works for VS14 and higher)
                cp -f ${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/parameters.log parameters_OLD.log
                
                runHubbing
                cat ${genomeName}/hub_address.txt
                echo
                
                rm -f ./${genomeName}/*.bw
                rm -f ./${genomeName}/*_statistics.log
                rm -rf ./READ*_fastqc_*
                rm -f ./read_trimming.log
                
                
            else
                # NOT GENERATING THIS AS IT IS DONE ALREADY IN THE PREVIOUS RUN AND IS WAITING IN THE PUBLIC FOLDER
                #ln -fs "${fileList[$i]}/${genomeName}/BigWigs/filtered_pileup_scaled.bw" .
                #echo "Generated symbolic link :"
                #echo "${fileList[$i]}/${genomeName}/BigWigs/filtered_pileup_scaled.bw"
                
                cd ..
                
                # For hubbing - info file copying.. (works for VS14 and higher)
                cp -f ${fileList[$i]}/${genomeName}/ForPeakCallAndFootPrintRerun/parameters.log parameters_OLD.log
                
                # Here take into account that for pyramid runs we only have one genome, and we don't have a folder for it..
                if [ "${PYRAMIDRERUN}" -eq "1" ]; then
                     genomeName="${TRUEgenomeName}"   
                fi
               
                runHubbing    
                
                cd  ${genomeName}
                cleanUpFolder
                listFolder
                echo
                cat hub_address.txt
                cd ..
                echo
                
            fi
            
            echo
            cd ..   
            
        done  



else

#--------THE-LOOP-over-all-ALREADY-MAPPED-bam-file-based-data-sets---------------------------------------

if [ "${SKIP_BOWTIE}" -eq 1 ] ; then
    
    folderList=($( cut -f 1 ./PIPE_mappedBamPaths.txt ))
    fileList=($( cut -f 2 ./PIPE_mappedBamPaths.txt ))
    GENOMEARRAY=($( cut -f 3 ./PIPE_mappedBamPaths.txt ))
    
    printRunStartArrays
    
#    cp -f ./PIPE_mappedBamPaths.txt ./dnase_pipe_1_param.txt
        
        #Run the rest of the pipe
        for i in $( seq 0 $((${#GENOMEARRAY[@]} - 1)) ); do
            genomeName=${GENOMEARRAY[$i]}
            printThis="Starting run for sample : ${folderList[$i]}, mapped to ${genomeName}"
            printToLogFile
            
            mkdir ${folderList[i]}
            mkdir ${folderList[i]}/${genomeName}
            
            cd ${folderList[i]}/${genomeName}
            pwd
            
            # Here have to circumvent a GRID ENGINE instability in ln command - sometimes LOSES underscores in this command and thus cannot complete what was asked.
            echo "ln -fs ${fileList[$i]} bowtie_out.bam" > temp.command
            chmod u=rwx temp.command
            ./temp.command
            rm -f temp.command
            echo "Generated symbolic link 'ln -fs ${fileList[$i]} bowtie_out.bam'"
            
            
            # Generate CO comment field in bam - this will hold the ORIGINAL READ ORDER which bowtie outputted (keeps read pairs with each others, if sorted along this flag)
            
            samtools view -o bowtie_out_noheading.sam bowtie_out.bam
            samtools view -H -o bowtie_out_heading.sam bowtie_out.bam
            
            
            # If we have comment field in the sam we read in - we delete it before appending a new one !
            
            coFieldPresent=-1
            coFieldPresent=$( head -n 1 bowtie_out_noheading.sam | awk '{for (a=1; a<=NF; a++) if (substr($a,1,5) == "CO:Z:") printf "%i" ,a ; else printf "%i" ,"-1" ;}' | sed 's/-1.*/-1/')
            
            if [ "${coFieldPresent}" -ne "-1" ]
            then
            cut -f ${coFieldPresent} --complement bowtie_out_noheading.sam > temp.sam
            mv -f temp.sam bowtie_out_noheading.sam
            fi
            
            # Combine sams to deed into comment field generation..
            cat bowtie_out_heading.sam bowtie_out_noheading.sam | grep -v "^$" > bowtie_out.sam
            
            testedFile="bowtie_out.sam"
            doTempFileTesting
            unlink bowtie_out.bam
            rm -f bowtie_out_heading.sam bowtie_out_noheading.sam
            
            generateCommentFieldAndSaveInBam
            
            
            cd ..
            
            setUCSCgenomeSizes
            
            runBowtieSubscript
            
            cd  ${genomeName}
            cleanAfterSam2bw
            cd ..
            
            runHubbing
            
            cd  ${genomeName}
            cleanUpFolder
            listFolder
            echo
            cat hub_address.txt
            cd ..   
            cd ..
        done    
else

#--------THE-LOOP-over-all-BAM-file-based-data-sets------------------------------------------------------

if [ -r "./PIPE_bamPaths.txt" ] ; then
    folderList=($( cut -f 1 ./PIPE_bamPaths.txt ))
    fileList=($( cut -f 2 ./PIPE_bamPaths.txt ))
    
    printRunStartArrays
    
    # This binary 1/0 tells if we originated to RunPipe from fastq or bam parameter file (and determines if the raw bam file is saved or not)
    fastqBasedRun=0    


for (( i=0; i<=$(( ${#folderList[@]} -1 )); i++ ))
do
    
    printThis="Starting run for sample : ${folderList[$i]}"
    printToLogFile
    
    mkdir ${folderList[$i]}
    cd ${folderList[$i]}
    pwd
    pwd >&2
    
    # Fetch BAM, transform to FASTQ :
    bamPath=${fileList[$i]}
    generateFastq
    # This is where BAM and FASTQ pipes merge
    runPipe
  
    cd ..

done

fi

#--------THE-LOOP-over-all-FASTQ-file-based-data-sets------------------------------------------------------



if [ -r "./PIPE_fastqPaths.txt" ] ; then
    
    
    folderList=($( cut -f 1 ./PIPE_fastqPaths.txt ))

    # Check how many columns we have.
    test=0
    if [ "${SINGLE_END}" -eq 0 ] ; then  
    test=$( cut -f 4 ./PIPE_fastqPaths.txt | grep -vc "^\s*$" )
    else
    test=$( cut -f 3 ./PIPE_fastqPaths.txt | grep -vc "^\s*$" )
    fi

    # If we have 3 columns paired end, or 2 columns single end :
    if [ "${test}" -eq "0" ]; then

    fileList1=($( cut -f 2 ./PIPE_fastqPaths.txt ))
    
    if [ "${SINGLE_END}" -eq 0 ] ; then
    fileList2=($( cut -f 3 ./PIPE_fastqPaths.txt ))
    fi
    
    # If we have 4 columns paired end, or 3 columns single end :
    else

    if [ "${SINGLE_END}" -eq 0 ] ; then
    cut -f 2,4 ./PIPE_fastqPaths.txt | awk '{ print $2"\t"$1 }' | tr "," "\t" | awk '{for (i=2;i<=NF;i++) printf "%s/%s,",$1,$i; print ""}' | sed 's/,$//' | sed 's/\/\//\//' > forRead1.txt
    cut -f 3,4 ./PIPE_fastqPaths.txt | awk '{ print $2"\t"$1 }' | tr "," "\t" | awk '{for (i=2;i<=NF;i++) printf "%s/%s,",$1,$i; print ""}' | sed 's/,$//' | sed 's/\/\//\//'  > forRead2.txt
    
    fileList1=($( cat ./forRead1.txt ))
    fileList2=($( cat ./forRead2.txt ))
    else
    cut -f 2,3 ./PIPE_fastqPaths.txt | awk '{ print $2"\t"$1 }' | tr "," "\t" | awk '{for (i=2;i<=NF;i++) printf "%s/%s,",$1,$i; print ""}' | sed 's/,$//' | sed 's/\/\//\//' > forRead1.txt
    fileList1=($( cat ./forRead1.txt ))
    fi
    
    rm -f forRead1.txt forRead2.txt
    
    fi
    
    printRunStartArraysFastq
    
    
    # This binary 1/0 tells if we originated to RunPipe from fastq or bam parameter file (and determines if the raw bam file is saved or not)
    fastqBasedRun=1


for (( i=0; i<=$(( ${#folderList[@]} -1 )); i++ ))
do
    printThis="Starting run for sample : ${folderList[$i]}"
    printToLogFile
    
    mkdir ${folderList[$i]}
    cd ${folderList[$i]}
    pwd
    pwd >&2
    
    # If we have single lane sequencing.
    if [ "$LANES" -eq 1 ] ; then 
    
    #Fetch FASTQ :
    fetchFastq
    
    else

    # If we have MULTIPLE lanes from sequencing.
    
    fetchFastqMultilane

    fi

    # This is where BAM and FASTQ pipes merge
    runPipe

    cd ..

done

#rm -f ./dnase_pipe_1b_param.txt



fi

#---------------------------------------

fi

fi

#--------Deleting-the-parameter-files-for-the-subscripts------------------------------------------------------

#rm -f parameters.log

if [ -r "./PIPE_hubbing.txt" ] ; then
   hubFolderCheck=$( cat ./PIPE_hubbing.txt | awk '{print $2"/"$1 }' )
elif [ -r "./PIPE_hubbingSymbolic.txt" ] ;  then
   hubFolderCheck=$( cat ./PIPE_hubbingSymbolic.txt | awk '{print $2"/"$1 }' )
fi

if [ -r "./PIPE_hubbing.txt" ] || [ -r "./PIPE_hubbingSymbolic.txt" ] ;  then 

#Remove the marker files from hub - run has ended !
rm -f ${hubFolderCheck}/HubFolder/*/${MagicNumber}.temp

rm -f ./dnase_pipe_2_param.txt
rm -f ./versionInfoHTML.txt

else
    
    echo
    echo "No hubbing done - as hubbing was not requested ! "

fi

timepoint=$( date )
echo "run finished : ${timepoint}"
exit 0


