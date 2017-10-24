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


# This is a trial to organize this script to 4 modules :
#
#########################################################
# MODULE 1) hub f1
#########################################################
# MODULE 2) hub f2 folder generation
#########################################################
# MODULE 3) hub f2 "other content" than bws - and their links to description.html
#########################################################
# MODULE 4) bigwigs - symbolic or not - and their links to tracks.txt 
#########################################################
# MODULE 5) listing the folders and hub address etc. 
#########################################################

doTrackExist(){
    # NEEDS THESE TO BE SET BEFORE CALL :
    #trackName=""
     
    if [ -s "${hubFolderHtsLongPath}/${genomeName}/tracks.txt" ]; then
    
    echo -e "grep -c \"track ${trackName}\$\" ${hubFolderHtsLongPath}/${genomeName}/tracks.txt" > temp.command
    chmod u=rwx temp.command
    trackExists=$(( $(./temp.command) ))
    rm -f temp.command
    
    else
    trackExists=0
    
    fi

}

writeBeginOfHtml(){
    
    echo "<!DOCTYPE HTML PUBLIC -//W3C//DTD HTML 4.01//EN" > begin.html
    echo "http://www.w3.org/TR/html4/strict.dtd" >> begin.html
    echo ">" >> begin.html
    echo " <html lang=en>" >> begin.html
    echo " <head>" >> begin.html
    echo " <title> ${hubName} data hub in ${genomeName} </title>" >> begin.html
    echo " </head>" >> begin.html
    echo " <body>" >> begin.html
    
}

doMultiWigParent(){
    
    # NEEDS THESE TO BE SET BEFORE CALL :
    #longLabel=""
    #trackName=""
    #overlayType=""
    #windowingFunction=""
    #visibility=""
    
    doTrackExist
    if [ "${trackExists}" -eq 0 ] ; then
    
    echo ""                                        >> ${genomeName}_tracks.txt
    echo "#--------------------------------------" >> ${genomeName}_tracks.txt
    echo ""                                        >> ${genomeName}_tracks.txt
    
    echo "track ${trackName}"                      >> ${genomeName}_tracks.txt
    echo "container multiWig"                      >> ${genomeName}_tracks.txt
    echo "shortLabel ${trackName}"                 >> ${genomeName}_tracks.txt
    echo "longLabel ${longLabel}"                  >> ${genomeName}_tracks.txt
    echo "type bigWig"                             >> ${genomeName}_tracks.txt
    echo "visibility ${visibility}"                >> ${genomeName}_tracks.txt
    #echo "aggregate transparentOverlay"           >> ${genomeName}_tracks.txt
    #echo "aggregate solidOverlay"                 >> ${genomeName}_tracks.txt
    echo "aggregate ${overlayType}"                >> ${genomeName}_tracks.txt
    echo "showSubtrackColorOnUi on"                >> ${genomeName}_tracks.txt
    #echo "windowingFunction maximum"              >> ${genomeName}_tracks.txt
    #echo "windowingFunction mean"                 >> ${genomeName}_tracks.txt
    echo "windowingFunction ${windowingFunction}"  >> ${genomeName}_tracks.txt
    echo "configurable on"                         >> ${genomeName}_tracks.txt
    echo "autoScale on"                            >> ${genomeName}_tracks.txt
    echo "alwaysZero on"                           >> ${genomeName}_tracks.txt
    echo "dragAndDrop subtracks"                   >> ${genomeName}_tracks.txt
    # Checking for double // in the path names - cases like this :  ${SERVERADDRESS}//public
    tempPath=$(echo "${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/description" | sed 's/\/\/*/\//g')
    echo "html ${SERVERTYPE}://${tempPath}"        >> ${genomeName}_tracks.txt
    echo ""                                        >> ${genomeName}_tracks.txt
    
    fi
}

doMultiWigChild(){
    
    # NEEDS THESE TO BE SET BEFORE CALL
    # parentTrack=""
    # trackName=""
    # fileName=".bw"
    # trackColor=""
    # trackPriority=""
    
    # Is this track already written to the tracks.txt file?
    doTrackExist
    if [ "${trackExists}" -eq 0 ] ; then
    
    # Does this track have data file which has non-zero size?
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/${fileName}" ]; then
    
    echo "track ${trackName}"                       >> ${genomeName}_tracks.txt
    echo "parent ${parentTrack}"                    >> ${genomeName}_tracks.txt
    # Checking for double // in the path names - cases like this :  ${SERVERADDRESS}//public
    tempPath2=$(echo "${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/${fileName}" | sed 's/\/\/*/\//g')
    echo "bigDataUrl ${SERVERTYPE}://${tempPath2}"  >> ${genomeName}_tracks.txt
    echo "shortLabel ${trackName}"                  >> ${genomeName}_tracks.txt
    echo "longLabel ${trackName}"                   >> ${genomeName}_tracks.txt
    echo "type bigWig"                              >> ${genomeName}_tracks.txt
    echo "color ${trackColor}"                      >> ${genomeName}_tracks.txt
    # Checking for double // in the path names - cases like this :  ${SERVERADDRESS}//public 
    tempPath=$(echo "${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/description" | sed 's/\/\/*/\//g')
    echo "html ${SERVERTYPE}://${tempPath}"         >> ${genomeName}_tracks.txt
    echo "priority ${trackPriority}"                >> ${genomeName}_tracks.txt
    echo ""                                         >> ${genomeName}_tracks.txt
    
    fi
    fi
}

doRegularTrack(){
    
    # NEEDS THESE TO BE SET BEFORE CALL
    # trackName=""
    # longLabel=""
    # fileName=".bw"
    # trackColor=""
    # trackPriority=""
    # visibility=""
    
    # Is this track already written to the tracks.txt file?
    doTrackExist
    if [ "${trackExists}" -eq 0 ] ; then
        
    # Does this track have data file which has non-zero size?
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/${fileName}" ] ; then

    echo ""                                          >> ${genomeName}_tracks.txt
    echo "#--------------------------------------"   >> ${genomeName}_tracks.txt
    echo ""                                          >> ${genomeName}_tracks.txt
    
    echo "track ${trackName}"                        >> ${genomeName}_tracks.txt
    # Checking for double // in the path names - cases like this :  ${SERVERADDRESS}//public
    tempPath2=$(echo "${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/${fileName}" | sed 's/\/\/*/\//g')
    echo "bigDataUrl ${SERVERTYPE}://${tempPath2}"   >> ${genomeName}_tracks.txt
    echo "shortLabel ${trackName}"                   >> ${genomeName}_tracks.txt
    echo "longLabel ${longLabel}"                    >> ${genomeName}_tracks.txt
    echo "autoScale on"                              >> ${genomeName}_tracks.txt
    echo "alwaysZero on"                             >> ${genomeName}_tracks.txt
    echo "type bigWig"                               >> ${genomeName}_tracks.txt
    echo "color ${trackColor}"                       >> ${genomeName}_tracks.txt
    echo "visibility ${visibility}"                  >> ${genomeName}_tracks.txt
    # Checking for double // in the path names - cases like this :  ${SERVERADDRESS}//public 
    tempPath=$(echo "${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/description" | sed 's/\/\/*/\//g')
    echo "html ${SERVERTYPE}://${tempPath}"          >> ${genomeName}_tracks.txt
    echo "priority ${trackPriority}"                 >> ${genomeName}_tracks.txt
    echo ""                                          >> ${genomeName}_tracks.txt
    
    fi
    fi
    
}

doPeakCallTrack(){
    
    # Like regular track - but draws only region from {0-1}, and in FULL visibility mode takes only 20pixels
    
    # NEEDS THESE TO BE SET BEFORE CALL
    # trackName=""
    # longLabel=""
    # fileName=".bw"
    # trackColor=""
    # trackPriority=""
    # visibility=""
    
    # Is this track already written to the tracks.txt file?
    doTrackExist
    if [ "${trackExists}" -eq 0 ] ; then
    
    # Does this track have data file which has non-zero size?
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/${fileName}" ]; then
        
    echo ""                                          >> ${genomeName}_tracks.txt
    echo "#--------------------------------------"   >> ${genomeName}_tracks.txt
    echo ""                                          >> ${genomeName}_tracks.txt
    
    echo "track ${trackName}"                        >> ${genomeName}_tracks.txt
    # Checking for double // in the path names - cases like this :  ${SERVERADDRESS}//public
    tempPath2=$(echo "${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/${fileName}" | sed 's/\/\/*/\//g')
    echo "bigDataUrl ${SERVERTYPE}://${tempPath2}"   >> ${genomeName}_tracks.txt
    echo "shortLabel ${trackName}"                   >> ${genomeName}_tracks.txt
    echo "longLabel ${longLabel}"                    >> ${genomeName}_tracks.txt
    #echo "autoScale on"                             >> ${genomeName}_tracks.txt
    #echo "alwaysZero on"                            >> ${genomeName}_tracks.txt
    echo "type bigWig"                               >> ${genomeName}_tracks.txt
    echo "color ${trackColor}"                       >> ${genomeName}_tracks.txt
    echo "visibility ${visibility}"                  >> ${genomeName}_tracks.txt
    echo "maxHeightPixels 10:20:30"                  >> ${genomeName}_tracks.txt
    echo "viewLimits 0:1"                            >> ${genomeName}_tracks.txt
    # Checking for double // in the path names - cases like this :  ${SERVERADDRESS}//public 
    tempPath=$(echo "${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/description" | sed 's/\/\/*/\//g')
    echo "html ${SERVERTYPE}://${tempPath}"          >> ${genomeName}_tracks.txt
    echo "priority ${trackPriority}"                 >> ${genomeName}_tracks.txt
    echo ""                                          >> ${genomeName}_tracks.txt
    
    fi
    fi
}

doControlTrack(){
    
    # Like regular track - but as PATH takes the FULL PATH
    
    # NEEDS THESE TO BE SET BEFORE CALL
    # trackName=""
    # longLabel=""
    # fileName="path/to/control/file.bw"
    # trackColor=""
    # trackPriority=""
    # visibility=""
    # windowingFunction=""
    
    # Is this track already written to the tracks.txt file?
    doTrackExist
    if [ "${trackExists}" -eq 0 ] ; then
    
    # Does this track have data file which has non-zero size?
    if [ -s "${fileName}" ] ; then

    echo ""                                          >> ${genomeName}_tracks.txt
    echo "#--------------------------------------"   >> ${genomeName}_tracks.txt
    echo ""                                          >> ${genomeName}_tracks.txt
    
    echo "track ${trackName}"                        >> ${genomeName}_tracks.txt
    # Checking for double // in the path names - cases like this :  ${SERVERADDRESS}//public
    tempPath2=$(echo "${SERVERADDRESS}/${fileName}" | sed 's/\/\/*/\//g')
    echo "bigDataUrl ${SERVERTYPE}://${tempPath2}"   >> ${genomeName}_tracks.txt
    echo "shortLabel ${trackName}"                   >> ${genomeName}_tracks.txt
    echo "longLabel ${longLabel}"                    >> ${genomeName}_tracks.txt
    
    echo "autoScale off"                             >> ${genomeName}_tracks.txt
    echo "viewLimits 0:9"                            >> ${genomeName}_tracks.txt
    echo "windowingFunction ${windowingFunction}"    >> ${genomeName}_tracks.txt
    
    echo "alwaysZero on"                             >> ${genomeName}_tracks.txt
    echo "type bigWig"                               >> ${genomeName}_tracks.txt
    echo "color ${trackColor}"                       >> ${genomeName}_tracks.txt
    echo "visibility ${visibility}"                  >> ${genomeName}_tracks.txt
    # Checking for double // in the path names - cases like this :  ${SERVERADDRESS}//public 
    tempPath=$(echo "${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/description" | sed 's/\/\/*/\//g')
    echo "html ${SERVERTYPE}://${tempPath}"          >> ${genomeName}_tracks.txt
    echo "priority ${trackPriority}"                 >> ${genomeName}_tracks.txt
    echo ""                                          >> ${genomeName}_tracks.txt
    
    fi
    fi
}

##################################################################################
#                                                                                #
# ABOVE = the FUNCTION DEFINITIONS (in bash they have to be first in the code)   #
# BELOW = the USAGE information (how to use the code)                            #
# IN THE END = the MAIN CODE itself.                                             #
#                                                                                #
##################################################################################


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

PipeTopPath="$( which $0 | sed 's/\/bin\/mainScripts\/dataHubGenerator.sh$//' )"

BashHelpersPath="${PipeTopPath}/bin/bashHelpers"

# PRINTING TO LOG AND ERROR FILES
. ${BashHelpersPath}/logFilePrinter.sh

#------------------------------------------

# From where to call the CONFIGURATION script..

confFolder="${PipeTopPath}/conf"

#------------------------------------------

echo
echo "PipeTopPath ${PipeTopPath}"
echo "confFolder ${confFolder}"
echo "BashHelpersPath ${BashHelpersPath}"
echo

#------------------------------------------

# Calling in the CONFIGURATION script and its default setup :


echo "Calling in the conf/serverAddressAndPublicDiskSetup.sh script and its default setup .."

SERVERTYPE="UNDEFINED"
SERVERADDRESS="UNDEFINED"
REMOVEfromPUBLICFILEPATH="NOTHING"
ADDtoPUBLICFILEPATH="NOTHING"
tobeREPLACEDinPUBLICFILEPATH="NOTHING"
REPLACEwithThisInPUBLICFILEPATH="NOTHING"

. ${confFolder}/serverAddressAndPublicDiskSetup.sh

setPublicLocations

echo
echo "SERVERTYPE ${SERVERTYPE}"
echo "SERVERADDRESS ${SERVERADDRESS}"
echo "ADDtoPUBLICFILEPATH ${ADDtoPUBLICFILEPATH}"
echo "REMOVEfromPUBLICFILEPATH ${REMOVEfromPUBLICFILEPATH}"
echo "tobeREPLACEDinPUBLICFILEPATH ${tobeREPLACEDinPUBLICFILEPATH}"
echo "REPLACEwithThisInPUBLICFILEPATH ${REPLACEwithThisInPUBLICFILEPATH}"
echo

#------------------------------------------


genomeName=""
magicNumber="RUN"
windowSize=3
pipePath="/t1-home/molhaem2/telenius/Jelena_DNase_pipe"
RERUN=0
SINGLEEND=0
PLOIDYFILTERED=1
SYMBOLICLINKS=0
ONLYHUB=0
OUTFILENAME="qsub.out"
ERRFILENAME="qsub.err"
ORANGEBLUE="0"

OPTS=`getopt -o o:,e:,g:,n:,w:,p:,W: --long rerun:,ploidyTrack:,onlyhub:,singleEnd:,symbolic:,orangeblue: -- "$@"`
if [ $? != 0 ]
then
    exit 1
fi

eval set -- "$OPTS"

while true ; do
    case "$1" in
        -g) genomeName=$2 ; shift 2;;
        -e) ERRFILENAME=$2 ; shift 2;;
        -o) OUTFILENAME=$2 ; shift 2;;
        -n) magicNumber=$2 ; shift 2;;
        -w) windowSize=$2 ; shift 2;;
        -W) windowTrackSize=$2 ; shift 2;;
        -p) pipePath=$2 ; shift 2;;
        --singleEnd) SINGLEEND=$2 ; shift 2;;
        --rerun) RERUN=$2 ; shift 2;;
        --onlyhub) ONLYHUB=$2 ; shift 2;;
        --ploidyTrack) PLOIDYFILTERED=$2 ; shift 2;;
        --symbolic) SYMBOLICLINKS=$2 ; shift 2;;
        --orangeblue) ORANGEBLUE=$2 ; shift 2;;
        --) shift; break;;
    esac
done

PloidyString="no ploidy regions, "
if [ "${PLOIDYFILTERED}" -eq 0 ]; then
    PloidyString=""
fi

echo "Starting run with parameters :"
echo "OUTFILENAME ${OUTFILENAME}"
echo "ERRFILENAME ${ERRFILENAME}"
echo "genomeName ${genomeName}"
echo "magicNumber ${magicNumber}"
echo "windowSize ${windowSize}"
echo "pipePath ${pipePath}"
echo "SINGLEEND ${SINGLEEND}"
echo "RERUN ${RERUN}"
echo "ONLYHUB ${ONLYHUB}"
echo "BLACKLISTFILTERED ${PLOIDYFILTERED}"
echo "SYMBOLICLINKS ${SYMBOLICLINKS}"
echo "ORANGEBLUE ${ORANGEBLUE}"
echo

# Setting track colors : normal or redgreen blindness ?

# RED GREEN (default)
redcolor="255,0,0"
greencolor="0,200,0"

# ORANGE BLUE (redgreen blindness colors)

if [ "${ORANGEBLUE}" -eq 1 ];then
    
redcolor="255,128,0"
greencolor="0,0,153"

fi


#########################################################
# MODULE 0.1) setting paths , setting timestamp
#########################################################

#-----------------------------------------------------------------

folderName=""

if [ -r "./dnase_pipe_1_param.txt" ] ; then
    
folderAsList=($( cut -f 1 ./dnase_pipe_1_param.txt ))
folderName=${folderAsList[0]}


else
    echo ""
    echo 'No parameter file dnase_pipe_1_param.txt provided in the run folder - pipeline aborted'
    echo ""
    exit 1
    
fi


#----------UCSC-data-generation-starts--------------------

if [ -r "../dnase_pipe_2_param.txt" ] ; then
    
    # Format of the dnase_pipe_2_param.txt
    #
    # MyHubName /public/telenius/MyNewDataHub  /public/telenius/MyfolderForBigWigs
    #
    
    # Generating lists for parameters. Of course we have only ONE parameter for both variables,
    # but we need to index it out with index [0] , as this construct generates list.

    hubNameList=($( cut -f 1 "../dnase_pipe_2_param.txt" ))
    hubName=${hubNameList[0]}
    
    hubFolderList=($( cut -f 2 "../dnase_pipe_2_param.txt" ))
    hubFolder=${hubFolderList[0]}
    
    # Here, parsing the data area location, to reach the public are address..
    
    diskFolder=${hubFolder}
    serverFolder=""
    
    echo
    parsePublicLocations
    echo
    
    # Setting the paths for the rest of the script ..
    
    hubFolderHtsLongPath="${diskFolder}/${hubName}/HubFolder"
    bigWigHtsLongPath="${diskFolder}/${hubName}/BigWigFolder"
    
    hubPublicPath="${serverFolder}/${hubName}/HubFolder"
    bigWigPublicShortPath="${serverFolder}/${hubName}/BigWigFolder"
    
    
    if [ "${SYMBOLICLINKS}" -eq 1 ] ; then
    
    # Save the symbolic link location
    symbolicStorageList=($( cut -f 3 "../dnase_pipe_2_param.txt" ))  
    realBigwigsForSymbolicPath="${symbolicStorageList[0]}/${hubName}"
    
    fi
    
    echo "hubPublicPath ${hubPublicPath}"
    echo "hubFolderHtsLongPath ${hubFolderHtsLongPath}"
    echo "bigWigPublicShortPath ${bigWigPublicShortPath}"
    echo "bigWigHtsLongPath ${bigWigHtsLongPath}"
    
    if [ "${SYMBOLICLINKS}" -eq 1 ] ; then
      
    echo "realBigwigsForSymbolicPath ${realBigwigsForSymbolicPath}"
        
    fi
    
#----- Generating TimeStamp - if files need backupping------------------

    TimeStamp=($( date | sed 's/[: ]/_/g' ))
    DateTime="$(date)"
    
    
#########################################################
# MODULE 0.2) checking file existence - and non-zero bigwig size 
#########################################################

#-----------------------------------------------------------------

# Checking for OBLIGATORY files and folders:

if [ "${RERUN}" -eq 0 ] ; then

#neededFiles[0]="./${genomeName}/mapped_statistics.log"
#neededFiles[1]="./${genomeName}/filtered_statistics.log"
#neededFiles[2]="./${genomeName}/bowtie_out_mapped_pileup.bw"
#neededFiles[3]="./${genomeName}/filtered_pileup.bw"

# neededFiles[0]="./${genomeName}/filtered_statistics.log"

neededFiles[0]="./${genomeName}/bowtie_out_mapped_pileup.bw"
neededFiles[1]="./${genomeName}/filtered_pileup.bw"


for (( f=0; f<=$(( ${#neededFiles[@]} -1 )); f++ ))
do
    
if [ ! -r "${neededFiles[$f]}" ] ; then
    echo "${neededFiles[$f]} file missing - ABORTING !" >&2
    exit 1
fi

if [ ! -s "${neededFiles[$f]}" ] ; then
    echo "${neededFiles[$f]} empty file - ABORTING !" >&2
    exit 1
fi
    
done

fi

#-----------------------------------------------------------------

    # Checking for files with zero size..
    
    for file in ./${genomeName}/*.bw
    do
        if [ ! -s "${file}" ]
        then
            rm -f "${file}"
        fi
    done
    
    #echo "Bigwig files with non-zero size.. :"
    #ls -l ./${genomeName}/*.bw


#########################################################
# MODULE 1) hub f1
#########################################################

    # If the Hub itself is not created already, we create it now..
    if [ ! -d "${hubFolderHtsLongPath}" ]; then
        
      printThis="Generating HUB folder ${hubFolderHtsLongPath}"
      printToLogFile   
      
      mkdir -p "${hubFolderHtsLongPath}"
      
      echo "hub ${hubName}" > TEMP_hub.txt
      echo "shortLabel ${hubName}" >> TEMP_hub.txt
      echo "longLabel ${hubName}" >> TEMP_hub.txt
      echo "genomesFile genomes.txt" >> TEMP_hub.txt
      echo "email jelena.telenius@gmail.com" >> TEMP_hub.txt
      mv TEMP_hub.txt "${hubFolderHtsLongPath}"/hub.txt
    fi
    
    if [ ! -d "${hubFolderHtsLongPath}/${genomeName}" ]; then
      
      mkdir "${hubFolderHtsLongPath}"/${genomeName}
      echo "genome ${genomeName}" > TEMP_genomes.txt
      echo "trackDb ${genomeName}/tracks.txt" >> TEMP_genomes.txt
      echo "" >> TEMP_genomes.txt
      cat TEMP_genomes.txt >> "${hubFolderHtsLongPath}"/genomes.txt
      rm -f TEMP_genomes.txt
      
    fi



#########################################################
# MODULE 2) hub f2 folder generation - and symbolic link folder generation (if requested)
#########################################################

   
    # If the public storage folder does not exist, creating it ..
    if [ ! -d "${bigWigHtsLongPath}" ] ; then
        printThis="Generating public storage folder ${bigWigHtsLongPath}"
        printToLogFile
        mkdir -p "${bigWigHtsLongPath}"
    fi
    
    if [ ! -d "${bigWigHtsLongPath}/${folderName}" ] ; then
        mkdir "${bigWigHtsLongPath}/${folderName}"
        printThis="Generating public storage folder ${bigWigHtsLongPath}/${folderName}"
        printToLogFile
    fi
    
    if [ ! -d "${bigWigHtsLongPath}/${folderName}/${genomeName}" ] ; then
        mkdir "${bigWigHtsLongPath}/${folderName}/${genomeName}"
        printThis="Generating public storage folder ${bigWigHtsLongPath}/${folderName}/${genomeName}"
        printToLogFile
    fi
    
    
    
    if [ "${SYMBOLICLINKS}" -eq 1 ] ; then
        
    # If the symbolic links folder does not exist, creating it ..
    if [ ! -d "${realBigwigsForSymbolicPath}" ] ; then
        printThis="Generating BigWig storage folder ${realBigwigsForSymbolicPath}/${hubName}"
        printToLogFile
        mkdir -p "${realBigwigsForSymbolicPath}"
    fi
    
    if [ ! -d "${realBigwigsForSymbolicPath}/${folderName}" ] ; then
        mkdir "${realBigwigsForSymbolicPath}/${folderName}"
        printThis="Generating BigWig storage folder ${realBigwigsForSymbolicPath}/${folderName}"
        printToLogFile
    fi
    
    if [ ! -d "${realBigwigsForSymbolicPath}/${folderName}/${genomeName}" ] ; then
        mkdir "${realBigwigsForSymbolicPath}/${folderName}/${genomeName}"
        printThis="Generating BigWig storage folder ${realBigwigsForSymbolicPath}/${folderName}/${genomeName}"
        printToLogFile
    fi

    fi

#########################################################
# MODULE 3) hub f2 "other content" than bws - and their links to description.html
#########################################################

    # HTML file generation - the BEGINNING OF THE FILE
        
    echo "Beginning to create html file ${bigWigHtsLongPath}/${folderName}/${genomeName}/description.html (on the side of copying files).."
    
    writeBeginOfHtml

    echo "" > temp_description.html
    
    if [ "${RERUN}" -eq 1 ] ; then
    echo "<hr />" >> temp_description.html
    fi
    
    echo "<p>Data produced ${DateTime} </p>" >> temp_description.html
    cat versionInfoHTML.txt >> temp_description.html
 
    if [ "${RERUN}" -eq 0 ] ; then
    echo "<p><h4>Data files produced by the script :</h4>" >> temp_description.html
    echo "<li><em>UNTRIMMED_bowtie_out.bam</em> = bowtie mapped data, all FASTQ reads (also unmapped), no trimming / flashing before mapping </li>" >> temp_description.html
    echo "<li><em>UNTRIMMED_onlyMapped_bowtie_out.bam</em> = bowtie mapped data, all FASTQ reads (only mapped), no trimming / flashing before mapping </li>" >> temp_description.html
    echo "<li><em>bowtie_out.bam</em> = bowtie mapped data (mapped after requested trimming / flashing) </li>" >> temp_description.html
    echo "<li><em>mapped.bam</em> = bowtie mapped data (mapped after requested trimming / flashing), only mapped reads (single end) / proper pairs (mapped end)</li>" >> temp_description.html
    echo "<li><em>filtered.bam</em> = bowtie mapped data, only mapped reads / proper pairs, ${PloidyString}no duplicates</li>" >> temp_description.html
    echo "<li><em>unpaired.bam</em> = bowtie mapped data, only mapped reads / proper pairs, ${PloidyString}no duplicates</li>" >> temp_description.html
    echo "<li><em>filtered_unpaired.bam</em> = bowtie mapped data, single end mapped after mapped end map failed, ${PloidyString}no duplicates</li>" >> temp_description.html
    echo "<li><em>LEFT.bdg</em> = footprinted mapped data, LEFTMOST 1base of the LEFT read</li>" >> temp_description.html
    echo "<li><em>RIGHT.bdg</em> = footprinted mapped data, RIGHTMOST 1base of the RIGHT read</li>" >> temp_description.html
    echo "</p>" >> temp_description.html
    echo "<hr />" >> temp_description.html
    fi
        
    printThis="Copying (non-bigwig) files to folder ${bigWigHtsLongPath} - OVERWRITES if existing files in storage location have same names"
    printToLogFile 
    
    if [ "${RERUN}" -eq 0 ] ; then
    
    rm -rf TEMPdirectory
    mkdir TEMPdirectory
    
    cp -fr READ*_fastqc_ORIGINAL TEMPdirectory/.  2> /dev/null
    cp -fr READ*_fastqc_PRE_TRIMMED TEMPdirectory/.  2> /dev/null
    cp -f ./read_trimming*.log TEMPdirectory/.  2> /dev/null
    
    mv -f TEMPdirectory/* ${bigWigHtsLongPath}/${folderName}/. 2> /dev/null
    
    rmdir TEMPdirectory
    
    rm -rf TEMPdirectory
    mkdir TEMPdirectory
    
    cp -fr ./${genomeName}/READ*_fastqc_TRIMMED TEMPdirectory/.  2> /dev/null
    
    cp -fr ./${genomeName}/UNMAPPED*_fastqc TEMPdirectory/.  2> /dev/null
    cp -fr ./${genomeName}/FLASHED*_fastqc TEMPdirectory/.  2> /dev/null
    
    cp -f ./${genomeName}/read_trimming*.log TEMPdirectory/.  2> /dev/null

    cp -f ./${genomeName}/flashing.log TEMPdirectory/.  2> /dev/null
    cp -f ./${genomeName}/flash.histogram TEMPdirectory/.  2> /dev/null
    cp -f ./${genomeName}/flash.hist TEMPdirectory/.  2> /dev/null

    cp -f ./${genomeName}/*_statistics.log TEMPdirectory/.  2> /dev/null
    cp -f ./${genomeName}/UNTRIMMED_bowtie_out_statistics.log TEMPdirectory/. 2> /dev/null
    cp -f ./${genomeName}/UNTRIMMED_onlyMapped_bowtie_out_statistics.log TEMPdirectory/.  2> /dev/null
    
    mv -f TEMPdirectory/* ${bigWigHtsLongPath}/${folderName}/${genomeName}/. 2> /dev/null
    
    rmdir TEMPdirectory 
    rm -rf TEMPdirectory   
    
    
    fi
    
    
    # Taking care of qsub file storage..
    
    originalQsubPath="EMPTY"
    rerunQsubPath="EMPTY"
    
    if [ "${RERUN}" -eq 0 ] && [ "${ONLYHUB}" -eq 0 ] ; then
    # If we are in "normal run" - we save "Original Files"
    
    cp -f ../${OUTFILENAME} ${bigWigHtsLongPath}/qsub_original.out 
    cp -f ../${ERRFILENAME} ${bigWigHtsLongPath}/qsub_original.err   
    
    
    hereIam=$( pwd ) 
    #Fetching the path..
    cd ..
    originalQsubPath=$( pwd )
    cd ${hereIam}

    else
    # We are a RERUN - so we need to check we have BOTH original AND rerun qsub files..
    
    # Trying to get original qsub files (if they are not saved here yet)
    
    # Trickstery to find the ORIGINAL running path (to print it out..)
    pathList=($( cut -f 2 ./dnase_pipe_1_param.txt ))
    # Saving where we are hubbing now..
    hereIam=$( pwd )
    #Fetching the path..
    cd ${pathList[0]}
    cd ..
    originalQsubPath=$( pwd )    
    # Going back..
    cd ${hereIam}

    # Going to fetch files..
    if [ ! -s "${bigWigHtsLongPath}/qsub_original.out" ]; then
        cp -f ${originalQsubPath}/${OUTFILENAME} ${bigWigHtsLongPath}/qsub_original.out 
    fi
    if [ ! -s "${bigWigHtsLongPath}/qsub_original.err" ]; then
        cp -f ${originalQsubPath}/${ERRFILENAME} ${bigWigHtsLongPath}/qsub_original.err 
    fi  
    
    
    # Saving the data of THIS RERUN, as well :
    # Generate "timestamp" and save to "all reruns qsub files" catenation..
    echo "" > echoDate
    echo "_RERUN_________________________________________" >> echoDate
    date >> echoDate
    echo "" >> echoDate
    
    cat echoDate ../${OUTFILENAME} > qsub_rerun_NEW.out
    cat echoDate ../${ERRFILENAME} > qsub_rerun_NEW.err
    
    # Catenate existing and new rerun qsub files (or generate from scratch)
    
    if [ -s "${bigWigHtsLongPath}/qsub_rerun.out" ] ; then
        
        cat qsub_rerun_NEW.out ${bigWigHtsLongPath}/qsub_rerun.out > temp.out
        mv -f temp.out ${bigWigHtsLongPath}/qsub_rerun.out
        
    else
        mv -f qsub_rerun_NEW.out ${bigWigHtsLongPath}/qsub_rerun.out
        
    fi
    
    if [ -s "${bigWigHtsLongPath}/qsub_rerun.err" ] ; then
        
        cat qsub_rerun_NEW.err ${bigWigHtsLongPath}/qsub_rerun.err  > temp.err
        mv -f temp.err ${bigWigHtsLongPath}/qsub_rerun.err
        
    else
        mv -f qsub_rerun_NEW.err ${bigWigHtsLongPath}/qsub_rerun.err
        
    fi
    
    
    # Fetching the path
    cd ..
    rerunQsubPath=$( pwd )    
    # Going back..
    cd ${hereIam}
    
    fi
    
    
    
    # HERE TO BE ADDED !!!
    
    # Here to be written the HTML file generation        

    echo "" >> temp_description.html    
    if [ "${RERUN}" -eq 0 ] ; then
    echo "<h3>Sample : ${folderName}</h3>" >> temp_description.html
    else
    echo "<h3>RERUN Sample : ${folderName}</h3>" >> temp_description.html
    fi    
    
    echo "<p>" >> temp_description.html
    #echo "Data located in : $(pwd)/${folderName}" >> temp_description.html
    echo "Data located in : $(pwd)" >> temp_description.html
    echo "</p>" >> temp_description.html
    echo "<p>" >> temp_description.html
    
   
    # If we have qsub files in ORIGINAL run
    if [ "${originalQsubPath}" != "EMPTY" ] ; then 
        # We assume that if the qsub.out file exists, also qsub.err exists
        if [ -s "${bigWigPublicShortPath}/qsub_original.out" ] ; then
          echo "<li>ORIGINAL RUN log files available in ${originalQsubPath}, as well as here : <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/qsub_original.out\" >${OUTFILENAME}</a> , and <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/qsub_original.err\" >${ERRFILENAME}</a>"  >> temp_description.html
        fi    
    fi
    
    # If we have qsub files in RERUN
    if [ "${rerunQsubPath}" != "EMPTY" ] ; then
        # We assume that if the qsub.out file exists, also qsub.err exists
        if [ -s ${bigWigPublicShortPath}/qsub_rerun.out ] ; then
          echo "<li>RERUN log files (this run) available in ${rerunQsubPath}, and combined log files of ALL the reruns done to this sample available here : <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/qsub_rerun.out\" >${OUTFILENAME}</a> , and <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/qsub_rerun.err\" >${ERRFILENAME}</a>"  >> temp_description.html
        fi
    fi
    
    echo "</p>" >> temp_description.html
    echo "<hr />" >> temp_description.html


    if [ "${RERUN}" -eq 0 ] ; then
    
    echo "<h4>Bowtie statistics here : </h4>" >> temp_description.html
    
    if [ -s ${bigWigPublicShortPath}/${folderName}/${genomeName}/UNTRIMMED_bowtie_out_statistics.log ] ; then
    echo "<li>\"FASTQ-to-bam\" bowtie run output - contains all reads from FASTQ : <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/UNTRIMMED_bowtie_out_statistics.log\" >UNTRIMMED_bowtie_out_statistics.log</a>  </li>" >> temp_description.html
    fi
    if [ -s ${bigWigPublicShortPath}/${folderName}/${genomeName}/UNTRIMMED_onlyMapped_bowtie_out_statistics.log ] ; then
    echo "<li>\"FASTQ-to-bam\" bowtie run output - contains all reads from FASTQ : <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/UNTRIMMED_onlyMapped_bowtie_out_statistics.log\" >UNTRIMMED_onlyMapped_bowtie_out_statistics.log</a>  </li>" >> temp_description.html
    fi
    echo "<li>Bowtie run output (unfiltered), all reads : <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/bowtie_out_statistics.log\" >bowtie_out_statistics.log</a>  </li>" >> temp_description.html
    echo "<li>Bowtie run output (unfiltered), only mapped reads : <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/mapped_statistics.log\" >mapped_statistics.log</a>  </li>" >> temp_description.html
    echo "<li>Bowtie run output, filtered (proper pairs, ${PloidyString}no duplicates) : <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/filtered_statistics.log\" >filtered_statistics.log</a>  </li>" >> temp_description.html
    
    UNPAIRED1=0
    UNPAIRED2=0
    UNPAIRED_FILTERED1=0
    UNPAIRED_FILTERED2=0

    # Checking if they have UNPAIRED reads..
    
    # UNPAIRED UNFILTERED-------------------------------------------------
    
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/singleEnd_READ1_statistics.log" ] || [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/singleEnd_READ2_statistics.log" ] ; then
    echo "<li>Unpaired reads (unfiltered) : " >> temp_description.html
    
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/singleEnd_READ1_statistics.log" ] ; then
        UNPAIRED1=1
        echo "<a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/singleEnd_READ1_statistics.log\" >READ1_statistics.log</a>  , and " >> temp_description.html
    fi
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/singleEnd_READ2_statistics.log" ] ; then
        UNPAIRED2=1
        echo "<a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/singleEnd_READ2_statistics.log\" >READ2_statistics.log</a>" >> temp_description.html
    fi
    
    echo "</li>" >> temp_description.html
    
    fi
    
    # UNPAIRED FILTERED-------------------------------------------------
    
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/singleEnd_READ1_filtered_statistics.log" ] || [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/singleEnd_READ2_filtered_statistics.log" ] ; then
    echo "<li>Unpaired filtered reads (${PloidyString}no duplicates) : " >> temp_description.html
    
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/singleEnd_READ1_filtered_statistics.log" ] ; then
        UNPAIRED_FILTERED1=1
        echo "<a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/singleEnd_READ1_filtered_statistics.log\" >READ1_filtered_statistics.log</a> , and ">> temp_description.html
    fi
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/singleEnd_READ2_filtered_statistics.log" ] ; then
        UNPAIRED_FILTERED2=1
        echo " <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/singleEnd_READ2_filtered_statistics.log\" >READ2_filtered_statistics.log</a>" >> temp_description.html
    fi
    
    echo "</li>" >> temp_description.html

    fi
    
    #--------------------------------------------------------------------

    # Checking if the fastqc was done for the samples :
    
    if [ -s "${bigWigHtsLongPath}/${folderName}/READ1_fastqc_ORIGINAL/fastqc_report.html" ] || [ -s "${bigWigHtsLongPath}/${folderName}/READ2_fastqc_ORIGINAL/fastqc_report.html" ]
    then
    echo "<h4>FASTQC results here : </h4>" >> temp_description.html

    echo "<li>FastQC results (untrimmed) : <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/READ1_fastqc_ORIGINAL/fastqc_report.html\" >READ1_fastqc_ORIGINAL/fastqc_report.html</a>   , and " >> temp_description.html
    echo " <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/READ2_fastqc_ORIGINAL/fastqc_report.html\" >READ2_fastqc_ORIGINAL/fastqc_report.html</a>  </li>" >> temp_description.html

    if [ -s "${bigWigHtsLongPath}/${folderName}/READ1_fastqc_PRE_TRIMMED/fastqc_report.html" ] || [ -s "${bigWigHtsLongPath}/${folderName}/READ2_fastqc_PRE_TRIMMED/fastqc_report.html" ]
    then
    TRIMMED=1
    echo "<li>FastQC results (pre-trim) - gently trimmed before bowtie : <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/READ1_fastqc_PRE_TRIMMED/fastqc_report.html\" >READ1_fastqc_PRE_TRIMMED/fastqc_report.html</a>  , and " >> temp_description.html
    echo " <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/READ2_fastqc_PRE_TRIMMED/fastqc_report.html\" >READ2_fastqc_PRE_TRIMMED/fastqc_report.html</a>  </li>" >> temp_description.html
    fi

    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/READ1_fastqc_TRIMMED/fastqc_report.html" ] || [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/READ2_fastqc_TRIMMED/fastqc_report.html" ]
    then
    TRIMMED=1
    echo "<li>FastQC results (final trim) - initially unmapped reads harshly trimmed and/or flashed :  <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/READ1_fastqc_TRIMMED/fastqc_report.html\" >READ1_fastqc_TRIMMED/fastqc_report.html</a>  , and " >> temp_description.html
    echo " <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/READ2_fastqc_TRIMMED/fastqc_report.html\" >READ2_fastqc_TRIMMED/fastqc_report.html</a>  </li>" >> temp_description.html
    fi
    
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/UNMAPPED_1_fastqc/fastqc_report.html" ] || [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/UNMAPPED_2_fastqc/fastqc_report.html" ]
    then
    echo "<li>FastQC results - only the reads not mapping (in proper pairs if PE data) even after trimming ( these reads were un-combinable via flash - if flash was requested) :  <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/UNMAPPED_1_fastqc/fastqc_report.html\" >UNMAPPED_1_fastqc/fastqc_report.html</a>  , and " >> temp_description.html
    echo " <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/UNMAPPED_2_fastqc/fastqc_report.html\" >UNMAPPED_2_fastqc/fastqc_report.html</a>  </li>" >> temp_description.html
    fi
    
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/FLASHED_1_fastqc/fastqc_report.html" ] 
    then
    echo "<li>FastQC results - the initially-unmapping, flash-combined reads :  <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/FLASHED_1_fastqc/fastqc_report.html\" >FLASHED_1_fastqc/fastqc_report.html</a> " >> temp_description.html
    fi
    
    fi
    
   
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/read_trimming_trim.log" ] || [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/flashing.log" ]
    then
    echo "<h4>Trimming/flashing log files here : </h4>" >> temp_description.html
    
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/read_trimming.log" ] 
    then
        echo "<li>Harsh trim_galore trim (for initially unmapped reads only) : <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/read_trimming.log\" >read_trimming.log</a>  </li>" >> temp_description.html
    fi
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/flashing.log" ]
    then
        echo "<li>Flashing (for initially unmapped reads only) : <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/flashing.log\" >flashing.log</a>  </li>" >> temp_description.html
    fi
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/flash.histogram" ]
    then
    echo "<li>Histogram of flashed reads : <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/flash.histogram\" >flash.histogram (\"picture\" format) </a> , <a target="_blank" href=\"${SERVERTYPE}://${SERVERADDRESS}/${bigWigPublicShortPath}/${folderName}/${genomeName}/flash.hist\" >flash.hist (\"spreadsheet\" format</a> </li>" >> temp_description.html
    fi
    
    fi
    
    echo "</p>"  >> temp_description.html
    echo "<hr />" >> temp_description.html
    

    # ending "if rerun=0 then write this part"    
    fi
    
    # The end of  HTML file 
    echo "</body>" > end.html
    echo "</html>"  >> end.html
    
    
    # Printing the files - depending on the situation (rerun, not rerun, overwrite, not overwrite etc)
    
    if [ "${RERUN}" -eq 0 ] ; then
    
    # NORMAL UPDATE - NORMAL RUN (NOT RERUN)
    
    # If the html file already exists, saving the old as copy (mv to a file name html_file_timestamp.html) 
    if [ -e "${bigWigHtsLongPath}/${folderName}/${genomeName}/description.html" ]
    then
        mv "${bigWigHtsLongPath}/${folderName}/${genomeName}/description.html" "${bigWigHtsLongPath}/${folderName}/${genomeName}/description.html_${TimeStamp}"
        echo "Generated backup for the stored html file : ${bigWigHtsLongPath}/${folderName}/${genomeName}/description.html_${TimeStamp}"
    fi
    
    cat begin.html temp_description.html end.html > "${bigWigHtsLongPath}/${folderName}/${genomeName}/description.html"
    rm -f begin.html temp_description.html end.html
    
    else

    # RERUN UPDATE - WE NEED TO CHECK IF WE WANT TO CREATE THE FILE FROM SCRATCH OR UPDATE EXISTING FILE

    if [ -e "${bigWigHtsLongPath}/${folderName}/${genomeName}/description.html" ]
    then
        cp "${bigWigHtsLongPath}/${folderName}/${genomeName}/description.html" "${bigWigHtsLongPath}/${folderName}/${genomeName}/description.html_${TimeStamp}"  2> /dev/null
        echo "Generated backup for the stored html file : ${bigWigHtsLongPath}/${folderName}/${genomeName}/description.html_${TimeStamp}"
        mv -f "${bigWigHtsLongPath}/${folderName}/${genomeName}/description.html" temp_file.html
        
        sed -i 's/<\/body>//' temp_file.html
        sed -i 's/<\/html>//' temp_file.html
        cat temp_file.html temp_description.html end.html > "${bigWigHtsLongPath}/${folderName}/${genomeName}/description.html"
    else
        cat begin.html temp_description.html end_description.html > "${bigWigHtsLongPath}/${folderName}/${genomeName}/description.html"
    fi
    
    rm -f begin.html temp_description.html end.html temp_file.html

    fi


    
#########################################################
# MODULE 4) bigwigs - symbolic or not - tracks.txt update - and their links to description.html
#########################################################  
    
    
    #----------------------------------------------------
    # Storing and linking bigwig files..
    
    if [ "${SYMBOLICLINKS}" -eq 1 ] ; then
    
    # Storing files in EVER-LASTING-NON-PUBLIC-COLDER
    
    printThis="Copying files to folder ${realBigwigsForSymbolicPath} - OVERWRITES if existing files in storage location have same names"
    printToLogFile
    rm -rf TEMPdirectory
    mkdir TEMPdirectory
    cp -f ./${genomeName}/*.bw TEMPdirectory/. 2> /dev/null
    mv -f TEMPdirectory/*.bw ${realBigwigsForSymbolicPath}/${folderName}/${genomeName}/.  2> /dev/null
    rmdir TEMPdirectory

    #echo "Bigwig files - now located in their permanent storage folder.. :"
    #ls -l ${realBigwigsForSymbolicPath}/${folderName}/${genomeName}/*.bw
    
    
    # Generating symbolic links
    
    printThis="Making symbolic links to folder ${bigWigHtsLongPath} - OVERWRITES if existing files in public folder have same names"
    printToLogFile
    
    thisDir=$( pwd )
    cd ${bigWigHtsLongPath}/${folderName}/${genomeName}
    ln -fs ${realBigwigsForSymbolicPath}/${folderName}/${genomeName}/*.bw .
    
    #echo "Bigwig files - now symbolically linked to their public folder.. :"
    #ls -l ${bigWigHtsLongPath}/${folderName}/${genomeName}/*.bw

    cd ${thisDir}
    
    else
    
    # Storing files in public folder.. 
    
    printThis="Copying files to folder ${bigWigHtsLongPath} - OVERWRITES if existing files in storage location have same names"
    printToLogFile 
    
    rm -rf TEMPdirectory
    mkdir TEMPdirectory
    cp -f ./${genomeName}/*.bw TEMPdirectory/. 2> /dev/null
    mv -f TEMPdirectory/*.bw ${bigWigHtsLongPath}/${folderName}/${genomeName}/. 2> /dev/null
    rmdir TEMPdirectory
    
    #echo "Bigwig files - now located in their public folder.. :"
    #ls -l ${bigWigHtsLongPath}/${folderName}/${genomeName}/*.bw
    
    fi
    
    
    
    #----------------------------------------------------
    # Generating the track files
    
    printThis="Generating tracks.txt file.."
    printToLogFile   
    
    # If the ${genomeName}_tracks.txt already exist, saving the old as copy (mv to a file name FileName_timestamp.txt)
    if [ -e "./${genomeName}_tracks.txt" ]
    then
        mv ${genomeName}_tracks.txt "${genomeName}_tracks.txt_${TimeStamp}" 2> /dev/null
    fi
    
    # Creating the output file with an empty line
    echo "" > ${genomeName}_tracks.txt
    
    #-----------------------------------------------------------------------------------------
        
    # Overlay track for the bowtie_out.bam and filtered.bam (to monitor how much got filtered)
    
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/bowtie_out_mapped_pileup.bw" ] || [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/filtered_pileup.bw" ]
    then
    
    # windowTrackSize
    longLabel="${folderName} all mapped reads RED, filtered reads GREEN"
    trackName="${folderName}"
    overlayType="solidOverlay"
    windowingFunction="maximum"
    visibility="full"
    doMultiWigParent

    # bowtie_out_mapped_pileup.bw
    # filtered_pileup.bw
    
    parentTrack="${folderName}"
    trackName="${folderName}_raw"
    fileName="bowtie_out_mapped_pileup.bw"
    trackColor=${redcolor}
    trackPriority="100"
    doMultiWigChild
    
    parentTrack="${folderName}"
    trackName="${folderName}_filtered"
    fileName="filtered_pileup.bw"
    trackColor=${greencolor}
    trackPriority="110"
    doMultiWigChild
    
    fi
    
    #-----------------------------------------------------------------------------------------
    
    # Overlay track WINDOWED for the bowtie_out.bam and filtered.bam (to monitor how much got filtered)
    
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/bowtie_out_mapped_window.bw" ] || [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/filtered_window.bw" ]
    then
    
    # windowTrackSize
    longLabel="${folderName} window size ${windowTrackSize}b, all mapped reads RED, filtered reads GREEN"
    trackName="${folderName}_window"
    overlayType="solidOverlay"
    windowingFunction="maximum"
    visibility="full"
    doMultiWigParent

    # bowtie_out_mapped_pileup.bw
    # filtered_pileup.bw
    
    parentTrack="${folderName}_window"
    trackName="${folderName}_window_raw"
    fileName="bowtie_out_mapped_window.bw"
    trackColor=${redcolor}
    trackPriority="100"
    doMultiWigChild
    
    parentTrack="${folderName}_window"
    trackName="${folderName}_window_filtered"
    fileName="filtered_window.bw"
    trackColor=${greencolor}
    trackPriority="110"
    doMultiWigChild
    
    fi

    # Overlay track for the UNPAIRED bowtie_out.bam and filtered.bam (to monitor how much got filtered)
    
    #--READ-1-----------------------------------------------------------------------------------------------
    if [ $((${UNPAIRED1})) -eq 1 ] || [ $((${UNPAIRED_FILTERED1})) -eq 1 ] ; then    
    longLabel="${folderName}, UNPAIRED read1 (mapped single end), - unfiltered RED, filtered (${PloidyString}no duplicates) GREEN"
    trackName="${folderName}_up1"
    overlayType="solidOverlay"
    windowingFunction="maximum"
    visibility="full"
    doMultiWigParent
       
    parentTrack="${folderName}_up1"
    trackName="${folderName}_raw_up1"
    fileName="singleEnd_READ1_pileup.bw" 
    trackColor=${redcolor}
    trackPriority="100"
    doMultiWigChild

    parentTrack="${folderName}_up1"
    trackName="${folderName}_up_filtered1"
    fileName="singleEnd_READ1_filtered_pileup.bw"
    trackColor=${greencolor}
    trackPriority="110"
    doMultiWigChild
    fi

    #--READ-2-----------------------------------------------------------------------------------------------
    if [ $((${UNPAIRED2})) -eq 1 ] || [ $((${UNPAIRED_FILTERED2})) -eq 1 ] ; then    
    longLabel="${folderName}, UNPAIRED read2 (mapped single end), - unfiltered RED, filtered (${PloidyString}no duplicates) GREEN"
    trackName="${folderName}_up2"
    overlayType="solidOverlay"
    windowingFunction="maximum"
    visibility="full"
    doMultiWigParent
       
    parentTrack="${folderName}_up2"
    trackName="${folderName}_raw_up2"
    fileName="singleEnd_READ2_pileup.bw" 
    trackColor=${redcolor}
    trackPriority="100"
    doMultiWigChild

    parentTrack="${folderName}_up2"
    trackName="${folderName}_up_filtered2"
    fileName="singleEnd_READ2_filtered_pileup.bw"
    trackColor=${greencolor}
    trackPriority="110"
    doMultiWigChild
    fi
   
    # Overlay track for the LEFT_pileup and RIGHT_pileup (i.e. the footprints)
    
    #-Normal-Bedgraphs-reads-------------------------------------------------------------------------------------
    
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/LEFT_pileup.bw" ] || [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/RIGHT_pileup.bw" ]
    then
    
    longLabel="${folderName} footprints (1st base), left-most read CYAN, right-most read ORANGE"
    trackName="${folderName}_footprint"
    overlayType="transparentOverlay"
    windowingFunction="maximum"
    visibility="full"
    doMultiWigParent
    
    # LEFT_pileup.bw
    # RIGHT_pileup.bw

    parentTrack="${folderName}_footprint"
    trackName="${folderName}_LEFT_foot"
    fileName="LEFT_pileup.bw"
    trackColor="90,170,255"
    trackPriority="100"
    doMultiWigChild

    parentTrack="${folderName}_footprint"
    trackName="${folderName}_RIGHT_foot"
    fileName="RIGHT_pileup.bw"
    trackColor="255,160,70"
    trackPriority="100"
    doMultiWigChild
    
    
    fi
    
    #-WINDOWED-reads-------------------------------------------------------------------------------------
    
    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/LEFT_window.bw" ] || [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/RIGHT_window.bw" ]
    then
        
    TEMPWINprint=$((${windowSize}*4))
    
    longLabel="${folderName} footprints (1st base), windowed, left CYAN, right ORANGE, readDepth (scaled by 0.1) GREY"
    trackName="${folderName}_footprintW"
    overlayType="transparentOverlay"
    windowingFunction="maximum"
    visibility="full"
    doMultiWigParent
    
    # LEFT_pileup.bw
    # RIGHT_pileup.bw
    
    parentTrack="${folderName}_footprintW"
    trackName="${folderName}_LEFT_footW"
    fileName="LEFT_window.bw"
    trackColor="90,170,255"
    trackPriority="100"
    doMultiWigChild
    
    parentTrack="${folderName}_footprintW"
    trackName="${folderName}_RIGHT_footW"
    fileName="RIGHT_window.bw"
    trackColor="255,160,70"
    trackPriority="100"
    doMultiWigChild

    parentTrack="${folderName}_footprintW"
    trackName="${folderName}_ReadsScaledby100"
    fileName="filtered_pileup_scaled.bw"
    trackColor="160,160,160"
    trackPriority="100"
    doMultiWigChild
    
    fi

    if [ -s "${bigWigHtsLongPath}/${folderName}/${genomeName}/LEFT_RIGHT_window.bw" ] 
    then
        
    longLabel="${folderName} footprints (1st base), windowed, left plus right Footprint BEIGE, readDepth (scaled by 0.01) GREY"
    trackName="${folderName}_footprintW2"
    overlayType="transparentOverlay"
    windowingFunction="maximum"
    visibility="full"
    doMultiWigParent
    
    parentTrack="${folderName}_footprintW2"
    trackName="${folderName}_LEFTRIGHT_footW"
    fileName="LEFT_RIGHT_window.bw"
    trackColor="220,200,160"
    trackPriority="100"
    doMultiWigChild


    parentTrack="${folderName}_footprintW2"
    trackName="${folderName}_ReadsScaledby100W2"
    fileName="filtered_pileup_scaled.bw"
    trackColor="160,160,160"
    trackPriority="100"
    doMultiWigChild
    
    fi
    
    
    #-PEAK-CALL-tracks----------------------------------------------------------------------------------------
    
       
    trackName="${folderName}_READ_peak"
    longLabel="${folderName} peak call - SeqMonk mimicking"
    fileName="READ_peakCall.bw"
    trackColor="0,0,200"
    trackPriority="100"
    visibility="full"
    
    #doRegularTrack
    doPeakCallTrack

  
    trackName="${folderName}_READ_peak_quant"
    longLabel="${folderName} peak call - SeqMonk mimicking"
    fileName="READ_peakCall_quant.bw"
    trackColor="0,0,200"
    trackPriority="100"
    visibility="hide"
    
    doRegularTrack
    #doPeakCallTrack

    
#---------------------------------------------------


    #-TRIMGALORE-control-track----------------------------------------------------------------------------------------
    
    #if [ $((${TRIMMED})) -eq 1 ] ; then
        
    #if [ ! -r  "/public/telenius/TRIMGALORE_for_pipe_290713/${genomeName}.bw" ] ; then
    #        echo "TrimmedAdapter track not available for the genome ${genomeName}" >&2
    #else 
    
    #trackName="TRIMGALORE_adapter"
    #longLabel="TrimGalore default adapter (from 4b AGAT to full lenght 13b AGATCGGAAGAGC) mapped in bowtie"
    #fileName="/public/telenius/TRIMGALORE_for_pipe_290713/${genomeName}.bw"
    #trackColor="0,0,0"
    #trackPriority="500"
    #visibility="dense"
    #windowingFunction="maximum"
    
    #doControlTrack
    
    #fi
    #fi
    
    #--------------------------------------------------------------------------------------------------------------------
    
    # If the ${genomeName}_tracks.txt already exist IN THE HUB, saving the old as copy (mv to a file name FileName_fimestamp.html)
    # This is naturally only done if the file existed before launching this pipeline (remember that multiple samples can update same tracks file!)
    
    if [ -e "${hubFolderHtsLongPath}/${genomeName}/tracks.txt" ]
    then
        
        # If we are not running yet (i.e. first time in this run we enter this directory) - we take safety copy - here magic number $$ is the shell id.
        if [ ! -r "${hubFolderHtsLongPath}/${genomeName}/${magicNumber}.temp" ]
          then
          cp "${hubFolderHtsLongPath}/${genomeName}/tracks.txt" "${hubFolderHtsLongPath}/${genomeName}/tracks.txt_${TimeStamp}"
        fi
        
        # In any case - we catenate what we have, to the existing tracks.
        mv -f "${hubFolderHtsLongPath}/${genomeName}/tracks.txt" previous_tracks.txt 2> /dev/null
        cat previous_tracks.txt ${genomeName}_tracks.txt > "${hubFolderHtsLongPath}/${genomeName}/tracks.txt"
        rm -f ${genomeName}_tracks.txt previous_tracks.txt
    
    # If the target directory had no tracks.txt
    else
        mv ${genomeName}_tracks.txt "${hubFolderHtsLongPath}/${genomeName}/tracks.txt"
    fi

    # Setting the magic marker - after first safety copy there should be THIS FILE in the folder - taking safety copies during run is prevented !
    echo "running!" > ${hubFolderHtsLongPath}/${genomeName}/${magicNumber}.temp
    # This file ${magicNumber}.temp is explicitly deleted in DnaseAndChip_pipe1.sh parent script, after the whole pipe has ran and finished.

    
#########################################################
# MODULE 5) listing the folders and hub address etc. 
#########################################################
    
   
    # Listing the produced files and data structures. Also giving the URL of the hub, and link to the "how to use a hub" pdf.
    
    echo
    echo "Data storage folder contents"
    echo "${bigWigHtsLongPath}/${folderName}"
    ls -lh "${bigWigHtsLongPath}/${folderName}" | cut -d "" -f 1,2,3,4 --complement 
    echo ""
    
    echo "Data storage folder contents - ${genomeName} folder"
    echo "${bigWigHtsLongPath}/${folderName}/${genomeName}"
    ls -lh "${bigWigHtsLongPath}/${folderName}/${genomeName}" | cut -d "" -f 1,2,3,4 --complement 
    echo ""
    echo "Generated / updated data hub :" > ${genomeName}/hub_address.txt
    
    tempPath=$(echo "${SERVERADDRESS}/${hubPublicPath}/hub.txt" | sed 's/\/\/*/\//g')
    
    echo "${SERVERTYPE}://${tempPath}" >> ${genomeName}/hub_address.txt
    echo >> ${genomeName}/hub_address.txt
   
    echo 'How to load this hub to UCSC : http://userweb.molbiol.ox.ac.uk/public/telenius/DataHubs/ReadMe/HowToUseA_DataHUB_160813.pdf' >> ${genomeName}/hub_address.txt
    cat ${genomeName}/hub_address.txt
    echo 'Other tutorials : http://userweb.molbiol.ox.ac.uk/public/telenius/DataHubs/ReadmeFiles.html' >> ${genomeName}/hub_address.txt

else
    printThis="No parameter file dnase_pipe_2_param.txt provided in the running folder.\n--> no data hub structures produced !"
    printToLogFile   
    
fi
