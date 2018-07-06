#!/bin/bash

##########################################################################
# Copyright 2017, Jelena Telenius (jelena.telenius@imm.ox.ac.uk)         #
#                                                                        #
# This file is part of CCseqBasic5 .                                     #
#                                                                        #
# CCseqBasic5 is free software: you can redistribute it and/or modify    #
# it under the terms of the GNU General Public License as published by   #
# the Free Software Foundation, either version 3 of the License, or      #
# (at your option) any later version.                                    #
#                                                                        #
# CCseqBasic5 is distributed in the hope that it will be useful,         #
# but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
# GNU General Public License for more details.                           #
#                                                                        #
# You should have received a copy of the GNU General Public License      #
# along with CCseqBasic5.  If not, see <http://www.gnu.org/licenses/>.   #
##########################################################################

#------------------------------------------

# Loading subroutines in ..

echo "Loading subroutines in .."

# /bin/runscripts/QC_and_Trimming.sh

PipeTopPath="$( echo $0 | sed 's/\/runscripts\/QC_and_Trimming.sh$//' )"

CapturePipePath="${PipeTopPath}/subroutines"

# DEBUG SUBROUTINES - for the situations all hell breaks loose
# . ${CapturePipePath}/debugHelpers.sh

# TESTING file existence, log file output general messages
CaptureCommonHelpersPath=$( dirname ${PipeTopPath} )"/commonSubroutines"
. ${CaptureCommonHelpersPath}/testers_and_loggers.sh
if [ "$?" -ne 0 ]; then
    printThis="testers_and_loggers.sh safety routines cannot be found in $0. Cannot continue without safety features turned on ! \n EXITING !! "
    printToLogFile
    exit 1
fi

#------------------------------------------

# From where to call the main scripts operating from the runscripts folder..

echo
echo "PipeTopPath ${PipeTopPath}"
echo "CapturePipePath ${CapturePipePath}"
echo


#------------------------------------------
# Setting $HOME to the current dir
echo 'Turning on safety measures for cd rm and mv commands in $0 - restricting script to file operations "from this dir below" only :'
HOME=$(pwd)
echo $HOME
echo
# - to enable us to track mv rm and cd commands
# from taking over the world accidentally
# this is done in testers_and_loggers.sh subroutines, which are to be used
# every time something gets parsed, and after that when the parsed value is used to mv cd or rm anything
# ------------------------------------------

# Test the testers and loggers ..

printThis="Testing the tester subroutines in $0 .."
printToLogFile
printThis="${CaptureCommonHelpersPath}/testers_and_loggers_test.sh 1> testers_and_loggers_test.out 2> testers_and_loggers_test.err"
printToLogFile
   
${CaptureCommonHelpersPath}/testers_and_loggers_test.sh 1> testers_and_loggers_test.out 2> testers_and_loggers_test.err
# The above exits if any of the tests don't work properly.

# The below exits if the logger test sub wasn't found (path above wrong or file not found)
if [ "$?" -ne 0 ]; then
    printThis="Testing testers_and_loggers.sh safety routines failed in $0. Cannot continue without testing safety features ! \n EXITING !! "
    printToLogFile
    exit 1
else
    printThis="Testing the tester subroutines completed - continuing ! "
    printToLogFile
fi

# Comment this out, if you want to save these files :
rm -f testers_and_loggers_test.out testers_and_loggers_test.err

#------------------------------------------

printThis="$0"
printToLogFile

# DOES NOT LOAD MODULES IT NEEDS :
#module load fastqc/0.11.4 
#module load trim_galore/0.3.1

# SHELL SCRIPTS inherit from their parents - (uncommenting) the following line proves that, if you are in doubt !
# module list 2>&1

# FILES IT REQUIRES :
# ${READ1}.fastq and ${READ2}.fastq - in the folder the script is launched.

# FILES IT OVERWRITES (if adapters were found) :
# ${READ1}.fastq and ${READ2}.fastq in the run directory.

# FILES AND FOLDERS IT PRODUCES (in the run directory) :

#read_trimming.log (if trimming needed)

# ${READ1}_fastqc_ORIGINAL
# ${READ2}_fastqc_ORIGINAL

# ${READ1}_fastqc_ORIGINAL.zip
# ${READ2}_fastqc_ORIGINAL.zip

# ${READ1}_fastqc_TRIMMED (if trimming needed)
# ${READ2}_fastqc_TRIMMED (if trimming needed) 

# ${READ1}_fastqc_TRIMMED.zip (if trimming needed)
# ${READ2}_fastqc_TRIMMED.zip (if trimming needed)


#---------------------------------------------------------------------

#echo "Code by Jelena Telenius, 24/10/2013"


#################################
# !!!!! OVERWRITE WARNING !!!!!
#################################

# NOTE THAT THE SCRIPT OVERWRITES THE FASTQ-FILES AS IT GOES 
# (ASSUMES THAT THERE IS SAVED BAM ORIGINAL FILES SOMEWHERE)
# BAM FILES ARE SMALLER THAN FASTQ FILES, SO IT IS BETTER TO SAVE THE READS LIKE THAT.
# IF YOU ONLY HAVE FASTQS, MAKE COPIES OF THEM BEFORE RUNNING THIS 
# (I.E. DO NOT RUN THIS TO YOUR ORIGINAL FILES)

# "Outline of the script :"
#
# Running fastqc
# "Outputting file encoding (illumina/sanger/?)"
# "Outputting overrepresented sequences info (pass/warn/fail) + the sequence and source"
# "If any adapter sequences (with suspected source) are found, trimming AND OVERWRITING the fastq files"
# "If trimming was needed, run fastqc again, to document that "now the data is trimmed better""
#


# READ NAMES are read in as flags $basename : used like this $basename.fastq
READ1="READ1"
READ2="READ2"
FLASHfastqNAME="flashed"
NEXTERA=0
CUSTOMAD=-1

# Default values Illumina PE adapters
A1_3prime="AGATCGGAAGAGC"
A2_3prime="AGATCGGAAGAGC"
A1_5prime="GCTCTTCCGATCT"
A2_5prime="GCTCTTCCGATCT"

replaceA1_3prime="no"
replaceA2_3prime="no"
replaceA1_5prime="no"
replaceA2_5prime="no"

# Quality score scheme for trim_galore..
QUAL=""

runMode=""
singleEnd=0

# trim_galore "defaults" for the pipeline :
stringency=1
length=10
qualFilter=20

OPTS=`getopt -o q: --long preFilter,filter5prime,flashFilter,fastqc,customad,qmin:,filter:,basenameR1:,basenameR2:,basename:,single:,nextera:,a31:,a32:,a51:,a52: -- "$@"`
if [ $? != 0 ]
then
    exit 1
fi

eval set -- "$OPTS"

while true ; do
    case "$1" in
        -q) QUAL=$2 ; shift 2;;
        --fastqc) runMode="fastqc" ; shift 1;;
        --customad) CUSTOMAD=1 ; shift 1;;
        --nextera) NEXTERA=$2 ; shift 2;;
        --qmin) qualFilter=$2 ; shift 2;;
        --filter) runMode=$2 ; shift 2;;
        --single) singleEnd=$2 ; shift 2;;
        --basenameR1) READ1=$2 ; shift 2;;
        --basenameR2) READ2=$2 ; shift 2;;
        --basename) FLASHfastqNAME=$2 ; shift 2;;
        --a31) replaceA1_3prime=$2 ; shift 2;;
        --a32) replaceA2_3prime=$2 ; shift 2;;
        --a51) replaceA1_5prime=$2 ; shift 2;;
        --a52) replaceA2_5prime=$2 ; shift 2;;
        --) shift; break;;
    esac
done

# NEXTERA adapters..
if [ "${NEXTERA}" -eq 1 ] ; then
    # There is one typo in this one (which was in use until 30 11 2015)
    #A1_3prime="CTGTCTCTTTT"
    #A2_3prime="CTGTCTCTTTT"
    #A1_5prime="AAAAGAGACAG"
    #A2_5prime="AAAAGAGACAG"
    
    # Corrected - to match the trim_galore new Nextera default (available from trim_galore VS 0.4.0 onwards)
    #Nextera:	CTGTCTCTTATA
    A1_3prime="CTGTCTCTTATA"
    A2_3prime="CTGTCTCTTATA"
    A1_5prime="TATAAGAGACAG"
    A2_5prime="TATAAGAGACAG"
fi

# Check if we override one or several adapters with custom ones..

if [ "${replaceA1_3prime}" != "no"  ] ; then
    A1_3prime=${replaceA1_3prime}
fi
if [ "${replaceA2_3prime}" != "no"  ] ; then
    A2_3prime=${replaceA2_3prime}
fi
if [ "${replaceA1_5prime}" != "no"  ] ; then
    A1_5prime=${replaceA1_5prime}
fi
if [ "${replaceA2_5prime}" != "no"  ] ; then
    A2_5prime=${replaceA2_5prime}
fi

echo "Starting run with parameters :" >&2
echo "runMode ${runMode}" >&2
echo "singleEnd ${singleEnd}" >&2
echo "READ1 ${READ1}" >&2
echo "READ2 ${READ2}" >&2
echo "FLASHfastqNAME ${FLASHfastqNAME}" >&2

if [ "${runMode}" != "fastqc" ]; then

echo "NEXTERA ${NEXTERA} (nextera adaptors = 1, illumina adaptors = 0) - may have been overwritten by separately given custom adaptors" >&2
echo "QUAL ${QUAL} ( quality score scheme )" >&2
echo "stringency ${stringency}" >&2
echo "length ${length}" >&2
echo "qualFilter ${qualFilter} ( quality filter treshold )" >&2
echo "A1_3prime ${A1_3prime}" >&2
echo "A2_3prime ${A2_3prime}" >&2
echo "A1_5prime ${A1_5prime}" >&2
echo "A2_5prime ${A2_5prime}" >&2

fi


# The above listed RUN MODES are the only flags of the script - they determine which part of this code is excecuted.

# The structure of the script is thus :

# if [ "${runMode}" = "filter" ]; then
#   RUN GREEDY FILTERING TRIM_GALORE (TO BE USED FOR NON-MAPPED BOWTIE READS)
#   exit 1
# fi
#
# if [ "${runMode}" = "fastq" ]; then
#   RUN ONLY FASTQC, NO FILTERING
#   exit 1
# fi
#

    

# check the input file existence.. (normal runs - 2 input files)
    if [ ! -r "./${READ1}.fastq" ] ; then
        printThis="Can not find file ${READ1}.fastq - aborting trimming"
        printToLogFile
        exit 1
    fi
    if [ ! -r "./${READ2}.fastq" ] && [ "${singleEnd}" -eq 0 ] ; then
        printThis="Can not find file ${READ2}.fastq - aborting trimming"
        printToLogFile
        exit 1
    fi


###################################################################
if [ "${runMode}" = "fastqc" ]; then
###################################################################
#   RUN ONLY FASTQC, NO FILTERING
###################################################################

echo "RUNNING FASTQC .."

echo "fastqc --quiet -f fastq ${READ1}.fastq"
fastqc -f fastq "${READ1}.fastq"

if [ "${singleEnd}" -eq 0 ] ; then
echo "fastqc --quiet -f fastq ${READ2}.fastq"
fastqc -f fastq "${READ2}.fastq"
fi

###################################################################
elif [ "${runMode}" -eq 3 ] ; then
###################################################################
#   RUN TRIM_GALORE - BASED ON ABOVE SET STRINGENCY AND LENGHT PARAMETERS
###################################################################
   
   S=${stringency}
   L=${length}
   A1=${A1_3prime}
   A2=${A2_3prime}
   
    printThis="Running trim_galore for ILLUMINA PAIRED END SEQUENCING ADAPTERS only - like this :\nTrims both reads against Illumina next-to-the-read 13bp adaptor sequence AGATCGGAAGAGC, (and everything in the read after it in 3' end)"
    printToTrimmingLogFile
   
    printThis="Allows 10% sequencing errors, at least ${S}b has to match, excludes both reads in pair if EITHER of them end up shorter than ${L}bp\nCuts the low-quality part of 3' end - quality cutoff being : ${qualFilter}, in PHRED${QUAL} quality score scheme"
    printToTrimmingLogFile

trimmingOK=1
if [ "${singleEnd}" -eq 0 ] ; then
    printThis="trim_galore --trim1 --phred${QUAL} --paired -q ${qualFilter} -a ${A1} -a2 ${A2} --length ${L} --stringency ${S} ${READ1}.fastq ${READ2}.fastq"
    printToTrimmingLogFile
    trim_galore --trim1 "--phred${QUAL}" --paired -q "${qualFilter}" -a "${A1}" -a2 "${A2}" --length "${L}" --stringency "${S}" "${READ1}.fastq" "${READ2}.fastq" >> "read_trimming.log"
    if [ $? -ne 0 ]; then trimmingOK=0;fi
else
    printThis="trim_galore --phred${QUAL} -q ${qualFilter} -a ${A1} --length ${L} --stringency ${S} ${READ1}.fastq"
    printToTrimmingLogFile
    trim_galore "--phred${QUAL}" -q "${qualFilter}" -a "${A1}" --length "${L}" --stringency "${S}" "${READ1}.fastq" >> "read_trimming.log"    
    if [ $? -ne 0 ]; then trimmingOK=0;fi
fi


    echo "" >> "read_trimming.log"
    echo "${READ1} - detailed trimming report :" >> "read_trimming.log"
    cat "${READ1}"*trimming_report.txt >> "read_trimming.log"
    
if [ "${singleEnd}" -eq 0 ] ; then
    echo "" >> "read_trimming.log"
    echo "${READ2} - detailed trimming report :" >> "read_trimming.log"
    cat "${READ2}"*trimming_report.txt >> "read_trimming.log"
    echo "" >> "read_trimming.log"
fi

    rm -f *trimming_report.txt
    
# Check if trimming went find, if not, exit 1

if [ "${trimmingOK}" -eq 0 ]; then
    printThis="Trimming crashed !\n EXITING !! "
    printToLogFile
    exit 1    
fi
    
# Check if files exist, if not, exit 1

if [ ! -s ${READ1}_val_1.fq ]; then
    printThis="Trimming failed : Trimgalore output file ${READ1}_val_1.fq is empty ! \n EXITING !! "
    printToLogFile
    exit 1    
fi

if [ "${singleEnd}" -eq 0 ] ; then
if [ ! -s ${READ2}_val_2.fq ]; then
    printThis="Trimming failed : Trimgalore output file ${READ2}_val_2.fq is empty ! \n EXITING !! "
    printToLogFile
    exit 1    
fi    
fi
    
    # OVERWRITE STEP !!!!
    
if [ "${singleEnd}" -eq 0 ] ; then
    
    moveCommand='mv -f ${READ1}_val_1.fq ${READ1}.fastq'
    moveThis="${READ1}"
    moveToHere="${READ1}"
    checkMoveSafety    
    mv -f "${READ1}_val_1.fq" "${READ1}.fastq"

    moveCommand='mv -f ${READ2}_val_1.fq ${READ2}.fastq'
    moveThis="${READ2}"
    moveToHere="${READ2}"
    checkMoveSafety      
    mv -f "${READ2}_val_2.fq" "${READ2}.fastq"
else
    moveCommand='mv -f ${READ1}_trimmed.fq ${READ1}.fastq'
    moveThis="${READ1}"
    moveToHere="${READ1}"
    checkMoveSafety  
    mv -f "${READ1}_trimmed.fq" "${READ1}.fastq" 
fi

echo
echo "TRIMMING STEP LOG FILE PRODUCED / UPDATED :"
ls "read_trimming.log"
echo ""

###################################################################
elif [ "${runMode}" -eq 5 ] ; then
###################################################################
#   RUN VERY GREEDY FRONT ADAPTER FILTERING (TO BE USED AFTER FLASH, TO REMOVE FRONT ADAPTER FROM THE READS)
###################################################################

   S=${stringency}
   L=${length}
   #A1="AGATCGGAAGAGC"
   #A2="AGATCGGAAGAGC"
   
   # Reverse
   # CGAGAAGGCTAGA
   # Complement
   # GCTCTTCCGATCT
   
   A1=${A1_5prime}
   A2=${A2_5prime}
   
    trimmingOK=1
   
    printThis="Running cutadapt for ILLUMINA PAIRED END SEQUENCING ADAPTERS only - like this :\nTrims both reads against the REVERSE COMPLEMENT (GCTCTTCCGATCT) of Illumina next-to-the-read 13bp adaptor sequence AGATCGGAAGAGC, (and everything before it in 5' end)"
    printToTrimmingLogFile
    
    printThis="cutadapt -O ${S} --quality-base=${QUAL} -q ${qualFilter} -m 0 -g ${A1} "
    printToTrimmingLogFile
    
    printThis="READ 1"
    printToTrimmingLogFile
    cutadapt -O "${S}" "--quality-base=${QUAL}" -q "${qualFilter}" -m 0 -g "${A1}" "${READ1}.fastq" -o TEMP_R1_trimmed.fastq >> "read_trimming.log"
    if [ $? -ne 0 ]; then trimmingOK=0;fi

if [ "${singleEnd}" -eq 0 ] ; then
    printThis="READ 2"
    printToTrimmingLogFile
    cutadapt -O "${S}" "--quality-base=${QUAL}" -q "${qualFilter}" -m 0 -g "${A1}" "${READ2}.fastq" -o TEMP_R2_trimmed.fastq >> "read_trimming.log"
    printThis="Running trim_galore in 'non-trimming mode' to eliminate too short read pairs, excludes both reads in pair if EITHER of them end up shorter than ${L}bp"
    printToTrimmingLogFile   
 
    printThis="trim_galore --trim1 --phred${QUAL} --paired -q ${qualFilter} -a ${A1} --length ${L} --stringency 1000 TEMP_R1_trimmed.fastq TEMP_R2_trimmed.fastq"
    printToTrimmingLogFile
    
    trim_galore "--phred${QUAL}" --trim1 --paired -q "${qualFilter}" -a "${A1}" --length "${L}" --stringency 1000 TEMP_R1_trimmed.fastq TEMP_R2_trimmed.fastq >> "read_trimming.log"
    if [ $? -ne 0 ]; then trimmingOK=0;fi
    rm -f TEMP_R1_trimmed.fastq TEMP_R2_trimmed.fastq

    moveCommand='mv -f TEMP_R1_trimmed_val_1.fq ${READ1}.fastq'
    moveThis="${READ1}"
    moveToHere="${READ1}"
    checkMoveSafety  
    mv -f "TEMP_R1_trimmed_val_1.fq" "${READ1}.fastq"
    mv -f "TEMP_R2_trimmed_val_2.fq" "${READ2}.fastq"
else
    moveCommand='mv -f TEMP_R1_trimmed.fastq ${READ1}.fastq'
    moveThis="${READ1}"
    moveToHere="${READ1}"
    checkMoveSafety  
    mv -f "TEMP_R1_trimmed.fastq" "${READ1}.fastq"
fi

if [ "${trimmingOK}" -eq 0 ]; then
    printThis="Trimming crashed !\n EXITING !! "
    printToLogFile
    exit 1    
fi
    
echo
echo "TRIMMING STEP LOG FILE PRODUCED / UPDATED :"
ls "read_trimming.log"
echo ""

###################################################################
else
    echo "QC_and_Trimming.sh runMode '${runMode}' not found !! - EXITING"
    exit 1

###################################################################
# Ending "main switch" if clause (list-pre-fastq-pre-normal-flash)
 fi
