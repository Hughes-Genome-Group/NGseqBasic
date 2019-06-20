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


echoTrimInfo(){

    trimInfo="Initially unmapped data will be"
    
    if [ "${TRIM}" -eq 1 ]; then
        trimInfo="${trimInfo} trimmed with trim_galore"
        if [ "${FLASH}" -eq 1 ]; then
            trimInfo="${trimInfo} and"
        fi
    fi
    if [ "${FLASH}" -eq 1 ]; then
        trimInfo="${trimInfo} ran through Flash to better map overlapping short reads"
    fi
    
    printThis="Preparing data for 'real bowtie run' in ${genomeName} build"
    printNewChapterToLogFile
    
    printThis="${trimInfo}"
    printToLogFile
    
    if [ "${TRIM}" -eq 0 ] && [ "${FLASH}" -eq 1 ]; then
    echo "NOTE !! INTERPRET YOUR RESULTS WITH CAUTION !! :"
    echo "Combination --flash --noTrim is NOT RECOMMENDED (can result in combining reads on the sites of ADAPTERS instead of the reads themselves). "
    echo
    fi
    
}

readCounts(){
    
    echo
    echo "Read counts :"
    
    mapped1count=0
    mapped2count=0
    unmapped1count=0
    unmapped2count=0
    filtered1count=0
    filtered2count=0
    
    if [ -s "MAPPED_1.fastq" ] ; then
        mapped1count=$(( $( grep -c "" MAPPED_1.fastq )/4 ))
    fi
    if [ -s "MAPPED_2.fastq" ] ; then
        mapped2count=$(( $( grep -c "" MAPPED_2.fastq )/4 ))
    fi
    
    if [ -s "UNMAPPED_1.fastq" ] ; then
        unmapped1count=$(( $( grep -c "" UNMAPPED_1.fastq )/4 ))
    fi
    if [ -s "UNMAPPED_2.fastq" ] ; then
        unmapped2count=$(( $( grep -c "" UNMAPPED_2.fastq )/4 ))
    fi
    
    if [ -s "M_FILTERED_1.fastq" ] ; then
        filtered1count=$(( $( grep -c "" M_FILTERED_1.fastq )/4 ))
    fi
    if [ -s "M_FILTERED_2.fastq" ] ; then
        filtered2count=$(( $( grep -c "" M_FILTERED_2.fastq )/4 ))
    fi
    
    echo "MAPPED_1.fastq : ${mapped1count}"
    if [ "${SINGLE_END}" -eq 0 ] ; then
    echo "MAPPED_2.fastq : ${mapped2count}"
    fi
    echo "UNMAPPED_1.fastq : ${unmapped1count}"
    if [ "${SINGLE_END}" -eq 0 ] ; then
    echo "UNMAPPED_2.fastq : ${unmapped2count}"
    fi
    echo "M_FILTERED_1.fastq : ${filtered1count}"
    if [ "${SINGLE_END}" -eq 0 ] ; then
    echo "M_FILTERED_2.fastq : ${filtered2count}"
    fi
    echo
    
}

trimBowtie(){
    
    # If there are reads remaining ..
    
    if [ -e READ1.fastq ] && [ -s READ1.fastq ] ; then
    
    echo "Preliminary bowtie to separate mapping and unmapping reads" >&2
    echo  >&2
    
    if [ "${BOWTIE}" -eq 2 ] ; then
    
    if [ "${SINGLE_END}" -eq 0 ] ; then
    echo "bowtie2 -p 3 ${bowtieQuals} --no-discordant --no-mixed -k 1 --maxins ${MAXINS} -x ${bowtieGenomeBuild} ${bowtieReadList} --un-conc UNMAPPED.fastq --al-conc MAPPED.fastq > /dev/null"
    bowtie2 -p 3 ${bowtieQuals} --no-discordant --no-mixed -k 1 --maxins ${MAXINS} -x ${bowtieGenomeBuild} ${bowtieReadList} --un-conc UNMAPPED.fastq --al-conc MAPPED.fastq > /dev/null
    rm -f READ*.fastq
    else
    echo "bowtie2 -p 3 ${bowtieQuals} -k 1 --maxins ${MAXINS} -x ${bowtieGenomeBuild} ${bowtieReadList} --un UNMAPPED.fastq --al MAPPED.fastq > /dev/null"
    bowtie2 -p 3 ${bowtieQuals} -k 1 --maxins ${MAXINS} -x ${bowtieGenomeBuild} ${bowtieReadList} --un UNMAPPED.fastq --al MAPPED.fastq > /dev/null
    rm -f READ*.fastq     
    fi

    elif [ "${BOWTIE}" -eq 3 ] ; then
    
    if [ "${SINGLE_END}" -eq 0 ] ; then
    echo "hisat2 -p 3 ${bowtieQuals} --no-spliced-alignment --no-discordant --no-mixed -k 1 --maxins ${MAXINS} -x ${bowtieGenomeBuild} ${bowtieReadList} --un-conc UNMAPPED.fastq --al-conc MAPPED.fastq > /dev/null"
    hisat2 -p 3 ${bowtieQuals} --no-spliced-alignment --no-discordant --no-mixed -k 1 --maxins ${MAXINS} -x ${bowtieGenomeBuild} ${bowtieReadList} --un-conc UNMAPPED.fastq --al-conc MAPPED.fastq > /dev/null
    rm -f READ*.fastq
    else
    echo "hisat2 -p 3 ${bowtieQuals} -k 1 --no-spliced-alignment --maxins ${MAXINS} -x ${bowtieGenomeBuild} ${bowtieReadList} --un UNMAPPED.fastq --al MAPPED.fastq > /dev/null"
    hisat2 -p 3 ${bowtieQuals} -k 1 --no-spliced-alignment --maxins ${MAXINS} -x ${bowtieGenomeBuild} ${bowtieReadList} --un UNMAPPED.fastq --al MAPPED.fastq > /dev/null
    rm -f READ*.fastq     
    fi

    else
    
    echo "bowtie -p 3 --chunkmb ${BOWTIEMEMORY} ${bowtieQuals} ${mParameter} --best --strata --maxins ${MAXINS} --sam ${bowtieGenomeBuild} ${bowtieReadList} --al MAPPED.fastq --un UNMAPPED.fastq --max M_FILTERED.fastq > /dev/null"
    bowtie -p 3 --chunkmb "${BOWTIEMEMORY}" ${bowtieQuals} ${mParameter} --best --strata --maxins ${MAXINS} --sam ${bowtieGenomeBuild} ${bowtieReadList} --al MAPPED.fastq --un UNMAPPED.fastq --max M_FILTERED.fastq > /dev/null
    rm -f READ*.fastq
   
    fi  
    
    
    # Now, we have at least SOME of these files :
    # MAPPED_1.fastq MAPPED_2.fastq
    # UNMAPPED_1.fastq UNMAPPED_2.fastq
    # M_FILTERED_1.fastq M_FILTERED_2.fastq
    
    if [ "${SINGLE_END}" -eq 1 ] ; then
        
       if [ -s "MAPPED.fastq" ] ; then 
       mv -f MAPPED.fastq MAPPED_1.fastq
       fi
       if [ -s "UNMAPPED.fastq" ] ; then 
       mv -f UNMAPPED.fastq UNMAPPED_1.fastq
       fi
       if [ -s "M_FILTERED.fastq" ] ; then 
       mv -f M_FILTERED.fastq M_FILTERED_1.fastq
       fi
       
    fi
    
    if [ "${BOWTIE}" -eq 2 ] &&  [ "${SINGLE_END}" -eq 0 ]; then

       if [ -s "MAPPED.1.fastq" ] ; then 
       mv -f MAPPED.1.fastq MAPPED_1.fastq
       fi
       if [ -s "MAPPED.2.fastq" ] ; then 
       mv -f MAPPED.2.fastq MAPPED_2.fastq
       fi
       if [ -s "UNMAPPED.1.fastq" ] ; then        
       mv -f UNMAPPED.1.fastq UNMAPPED_1.fastq
       fi
       if [ -s "UNMAPPED.2.fastq" ] ; then 
       mv -f UNMAPPED.2.fastq UNMAPPED_2.fastq
       fi

       
    fi
    
    if [ "${BOWTIE}" -eq 3 ] &&  [ "${SINGLE_END}" -eq 0 ]; then

       if [ -s "MAPPED.1.fastq" ] ; then 
       mv -f MAPPED.1.fastq MAPPED_1.fastq
       fi
       if [ -s "MAPPED.2.fastq" ] ; then 
       mv -f MAPPED.2.fastq MAPPED_2.fastq
       fi
       if [ -s "UNMAPPED.1.fastq" ] ; then        
       mv -f UNMAPPED.1.fastq UNMAPPED_1.fastq
       fi
       if [ -s "UNMAPPED.2.fastq" ] ; then 
       mv -f UNMAPPED.2.fastq UNMAPPED_2.fastq
       fi

       
    fi
    
    echo
    echo "Generated files :"
    # mappED unmappED filterED
    ls -lht  | grep "ED_[12].fastq" 
    
    readCounts
    
    else
        
    printThis="No unmapping reads remaining - skipping new bowtie run !"
    printToLogFile
    
    fi
    
    
}

trimRound(){
    
    #Needs this to be set before call
    #trimDirection=3 or 5 (3' or 5' direction trimming)
    
    # We have custom adapters
    if [ "${CUSTOMAD}" -eq 1 ] ; then
        
    printThis="${PipePath}/QC_and_Trimming.sh -q ${intQuals} --filter ${trimDirection} --basenameR1 UNMAPPED_1 --basenameR2 UNMAPPED_2 --single ${SINGLE_END} --nextera ${NEXTERA} --customad --a31 ${ADA31} --a32 ${ADA32} --a51 ${ADA51} --a52 ${ADA52}"
    printToLogFile
    ${PipePath}/QC_and_Trimming.sh -q "${intQuals}" --filter "${trimDirection}" --basenameR1 "UNMAPPED_1" --basenameR2 "UNMAPPED_2" --single ${SINGLE_END} --nextera ${NEXTERA} --customad --a31 ${ADA31} --a32 ${ADA32} --a51 ${ADA51} --a52 ${ADA52}
    if [ "$?" -ne 0 ]; then
    printThis="TrimGalore run failed ! Possible reasons : \n 1) did you maybe use .gz packed files without adding --gz to the run parameters ? \n 2) did you try to run with corrupted input fastq files ? \n EXITING !! "
    printToLogFile
    exit 1
    fi
        
    else
    # We have standard adapters, one way or another.
    
    printThis="${PipePath}/QC_and_Trimming.sh -q ${intQuals} --filter ${trimDirection} --basenameR1 UNMAPPED_1 --basenameR2 UNMAPPED_2 --single ${SINGLE_END} --nextera ${NEXTERA}"
    printToLogFile
    ${PipePath}/QC_and_Trimming.sh -q "${intQuals}" --filter "${trimDirection}" --basenameR1 "UNMAPPED_1" --basenameR2 "UNMAPPED_2" --single ${SINGLE_END} --nextera ${NEXTERA}
    # This updated UNMAPPED_1.fastq UNMAPPED_2.fastq to be the TRIMMED files (overwrote them with the trimmed versions)
    if [ "$?" -ne 0 ]; then
    printThis="TrimGalore run failed ! Possible reasons : \n 1) did you maybe use .gz packed files without adding --gz to the run parameters ? \n 2) did you try to run with corrupted input fastq files ? \n EXITING !! "
    printToLogFile
    exit 1
    fi
    
    fi
    
    echo
    
    # Renaming files to not to overwrite..
    
    
    if [ -s "MAPPED_1.fastq" ] ; then  
    mv -f MAPPED_1.fastq "MAPPED_pre${trimDirection}end_1.fastq"
    fi
    if [ -s "M_FILTERED_1.fastq" ] ; then  
    mv -f M_FILTERED_1.fastq "M_FILTERED_pre${trimDirection}end_1.fastq"
    fi
    
    if [ "${SINGLE_END}" -eq 0 ] ; then

    if [ -s "MAPPED_2.fastq" ] ; then  
    mv -f MAPPED_2.fastq "MAPPED_pre${trimDirection}end_2.fastq"
    fi
    if [ -s "M_FILTERED_2.fastq" ] ; then  
    mv -f M_FILTERED_2.fastq "M_FILTERED_pre${trimDirection}end_2.fastq"
    fi

    fi  
    
    # Renaming to become input for next bowtie run..
    
    if [ -s "UNMAPPED_1.fastq" ] ; then    
    mv -f UNMAPPED_1.fastq READ1.fastq
    fi
    
    if [ "${SINGLE_END}" -eq 0 ] ; then
    if [ -s "UNMAPPED_2.fastq" ] ; then
    mv -f UNMAPPED_2.fastq READ2.fastq
    fi
    fi

    # Bowtie round - what maps, maps, what not, will be marked "unmapped"
    printThis="Trying to map the PREVIOUSLY UNMAPPING reads after trimming.."
    printToLogFile
    
    trimBowtie
    
    doQuotaTesting
    
}

trimmingFlashing(){
    # This is the code which does ALL the trimming and flashing before "real bowtie run".
    # It needs READ1.fastq, READ2.fastq, and replaces those files with the trimmed/flashed versions.
    # It also runs fastqc after all the hassle, to provide proof that "changes were made".
    # THESE PARAMETERS NEED TO BE SET :
    # mParameter=
    # genomeName=
    
    
    # ROUND one - what maps, maps, what not, will be filtered.
    
    printThis="Separating unmapped reads from ALL READS.."
    printToLogFile
    
    trimBowtie
    
    # Now, we have at least SOME of these files :
    # MAPPED_1.fastq MAPPED_2.fastq
    # UNMAPPED_1.fastq UNMAPPED_2.fastq
    # M_FILTERED_1.fastq M_FILTERED_2.fastq
    
  
    # UNMAPPED_1.fastq UNMAPPED_2.fastq --> these will now continue to filtering and flashing (as was requested in the flags)
    
    if [ "${TRIM}" -eq 1 ] ; then
    
    echo "Trimming unmapped reads.."
    
    #################################
    # 3' end trimming
    #################################
    
    printThis="3' end trimming"
    printToLogFile
    
    trimHere=0
    
    # Check if trimming is needed - i.e. if there are unmapped reads..
    if [ "$SINGLE_END" -eq 1 ] ; then
     
        if [ ! -s UNMAPPED_1.fastq ] && [ "$SINGLE_END" -eq 1 ] ; then
        printThis="No unmapped reads found. Skipping trimming.."
        printToLogFile
        else
        trimHere=1
        fi
        
    else
    
        if [ ! -s UNMAPPED_1.fastq ] || [ ! -s UNMAPPED_2.fastq ] ; then   
        printThis="No unmapped read pairs found. Skipping trimming.."
        printToLogFile
        else
        trimHere=1
        fi         
    fi
    
    if [ "$trimHere" -eq 1 ] ; then   

    trimDirection=3
    trimRound
    
    fi

    #################################
    # 5' end trimming
    #################################
    
    if [ "$TRIM5" -eq 1 ] ; then
    
    printThis="5' end trimming"
    printToLogFile

    # Check if trimming is needed - i.e. if there are unmapped reads..
    if [ "$SINGLE_END" -eq 1 ] ; then
     
        if [ ! -s UNMAPPED_1.fastq ] && [ "$SINGLE_END" -eq 1 ] ; then
        printThis="No unmapped reads found. Skipping trimming.."
        printToLogFile
        else
        trimHere=1
        fi
        
    else
    
        if [ ! -s UNMAPPED_1.fastq ] || [ ! -s UNMAPPED_2.fastq ] ; then   
        printThis="No unmapped read pairs found. Skipping trimming.."
        printToLogFile
        else
        trimHere=1
        fi         
    fi
    
    if [ "$trimHere" -eq 1 ] ; then   

    trimDirection=5
    trimRound
    
    fi
    
    fi
    
    ######################################
    
    printThis="Combining mapped and m-filtered files.."
    printToLogFile
    
    # Now we have some (or all) of these files :
    
    #MAPPED_pre3end_1.fastq MAPPED_pre5end_1.fastq MAPPED_1.fastq 
    #MAPPED_pre3end_2.fastq MAPPED_pre5end_2.fastq MAPPED_2.fastq 
    #M_FILTERED_pre3end_1.fastq M_FILTERED_pre5end_1.fastq M_FILTERED_1.fastq 
    #M_FILTERED_pre3end_1.fastq M_FILTERED_pre5end_1.fastq M_FILTERED_1.fastq 
    #UNMAPPED_1.fastq UNMAPPED_2.fastq
    
    # Making catenation series (the files have to go in right order, so using * when combining is dangerous)
    
    mapped3prime1=""
    mapped3prime2=""
    mapped5prime1=""
    mapped5prime2=""
    mapped1=""
    mapped2=""

    filtered3prime1=""
    filtered3prime2=""
    filtered5prime1=""
    filtered5prime2=""
    filtered1=""
    filtered2=""
    
    # Assuming that if we have R1, we have equally many in R2..
    
    # MAPPED reads..
    if [ -s MAPPED_pre3end_1.fastq ] ; then
        mapped3prime1="MAPPED_pre3end_1.fastq"
        mapped3prime2="MAPPED_pre3end_2.fastq"  
    fi
    
    if [ -s MAPPED_pre5end_1.fastq ] ; then
        mapped5prime1="MAPPED_pre5end_1.fastq"
        mapped5prime2="MAPPED_pre5end_2.fastq"  
    fi
    
    if [ -s MAPPED_1.fastq ] ; then
        mapped1="MAPPED_1.fastq"
        mapped2="MAPPED_2.fastq"  
    fi
    
    cat ${mapped3prime1} ${mapped5prime1} ${mapped1} > mapped_1.fastq
    testedFile="mapped_1.fastq"
    doTempFileFYI
    rm -f MAPPED_*1.fastq
    
    if [ "${SINGLE_END}" -eq 0 ] ; then 
    cat ${mapped3prime2} ${mapped5prime2} ${mapped2} > mapped_2.fastq
    testedFile="mapped_2.fastq"
    doTempFileFYI
    rm -f MAPPED_*2.fastq
    fi
    
    cat mapped_1.fastq | grep -v "^$" > MAPPED_1.fastq
    rm -f mapped_1.fastq
    
    if [ "${SINGLE_END}" -eq 0 ] ; then    
    cat mapped_2.fastq | grep -v "^$" > MAPPED_2.fastq
    rm -f mapped_2.fastq
    fi
    
    # FILTERED reads..
    if [ -s M_FILTERED_pre3end_1.fastq ] ; then
        filtered3prime1="M_FILTERED_pre3end_1.fastq"
        filtered3prime2="M_FILTERED_pre3end_2.fastq"  
    fi
    
    if [ -s M_FILTERED_pre5end_1.fastq ] ; then
        filtered5prime1="M_FILTERED_pre5end_1.fastq"
        filtered5prime2="M_FILTERED_pre5end_2.fastq"  
    fi
    
    if [ -s M_FILTERED_1.fastq ] ; then
        filtered1="M_FILTERED_1.fastq"
        filtered2="M_FILTERED_2.fastq"  
    fi
    
    cat ${filtered3prime1} ${filtered5prime1} ${filtered1} > filtered_1.fastq
    testedFile="filtered_1.fastq"
    doTempFileFYI
    rm -f M_FILTERED_*1.fastq
    
    if [ "${SINGLE_END}" -eq 0 ] ; then
    cat ${filtered3prime2} ${filtered5prime2} ${filtered2} > filtered_2.fastq
    testedFile="filtered_2.fastq"
    doTempFileFYI
    rm -f M_FILTERED_*2.fastq
    fi
    
    cat filtered_1.fastq | grep -v "^$" > M_FILTERED_1.fastq
    rm -f filtered_1.fastq

    if [ "${SINGLE_END}" -eq 0 ] ; then    
    cat filtered_2.fastq | grep -v "^$" > M_FILTERED_2.fastq
    rm -f filtered_2.fastq
    fi
    
    echo
    echo "Counting total read counts after trimming and re-mapping.."
    
    readCounts
    
    # We exit trimming with files :
    
    # Combined originally-mapping and after-trimming-mapping :
    # MAPPED_1.fastq MAPPED_2.fastq
    
    # Combined not-yet-mapping :
    # Filtering1.fastq Filtering2.fastq
    
    else
    printThis="No trimming of unmapped reads requested - skipping fastq trimming !"
    printToLogFile

    #------------------------------------------
    # end of "if [ "${TRIM}" -eq 1 ] ; then"
    fi
    #------------------------------------------
    

    if [ "$SINGLE_END" -eq 1 ] ; then
    printThis="Single-end run! Skipping flashing."
    printToLogFile
    else

    if [ "${FLASH}" -eq 1 ] ; then
    if [ ! -s UNMAPPED_1.fastq ] || [ ! -s UNMAPPED_2.fastq ] ; then
    
    printThis="No unmapped reads left. Skipping flashing."
    printToLogFile
    
    else

    echo
    echo "Running flash with parameters :"
    echo " m (minimum overlap) 9b"
    echo " x (sum-of-mismatches/overlap-lenght) = 1/8 = 0.125"
    echo " p phred score min (33 or 64)"
    echo
    printThis="flash -m 9 -x 0.125 -p ${intQuals} UNMAPPED_1.fastq UNMAPPED_2.fastq > flashing.log"
    printToLogFile
    
    flash -m 9 -x 0.125 -p "${intQuals}" UNMAPPED_1.fastq UNMAPPED_2.fastq > flashing.log
    
    ls | grep out*fastq
    
    # This outputs these files :
    # flashing.log  out.extendedFrags.fastq  out.hist  out.histogram  out.notCombined_1.fastq  out.notCombined_2.fastq
    rm -f UNMAPPED_1.fastq UNMAPPED_2.fastq
    
    if [ -s "out.extendedFrags.fastq" ]; then
    mv -f out.extendedFrags.fastq extendedFrags.fastq
    fi
    if [ -s "out.notCombined_1.fastq" ]; then
    mv -f out.notCombined_1.fastq UNMAPPED_1.fastq  
    fi
    if [ -s "out.notCombined_2.fastq" ]; then
    mv -f out.notCombined_2.fastq UNMAPPED_2.fastq
    fi
    
    mv -f out.hist flash.hist
    mv -f out.histogram flash.histogram
    

    echo "Read counts after flash :"
    
    flashedCount=0
    unmapped1count=0
    unmapped2count=0
    
    if [ -s "extendedFrags.fastq" ] ; then
        flashedCount=$(( $( grep -c "" extendedFrags.fastq )/4 ))
    fi
    if [ -s "UNMAPPED_1.fastq" ] ; then
        unmapped1count=$(( $( grep -c "" UNMAPPED_1.fastq )/4 ))
    fi
    if [ -s "UNMAPPED_2.fastq" ] ; then
        unmapped2count=$(( $( grep -c "" UNMAPPED_2.fastq )/4 ))
    fi
    
    echo "extendedFrags.fastq (count of read pairs combined in flash) : ${flashedCount}"
    echo "UNMAPPED_1.fastq (not extendable via flash) : ${unmapped1count}"
    echo "UNMAPPED_2.fastq (not extendable via flash) : ${unmapped2count}"
    
    
    if [ -s extendedFrags.fastq ] ; then
    
    printThis="Generating PE fragments from flash-combined reads.."
    printToLogFile
    
    # Running reverse complementing
    #- generating file where read names have been changed from R1 names to R2 names,
    # and Qual reversed, and Seq rev complemented :
    echo "perl ${PerlHelpersPath}/reverse_seq.pl extendedFrags.fastq" >&2
    perl ${PerlHelpersPath}/reverse_seq.pl extendedFrags.fastq
    # returned file extendedFrags_RevCompl.fastq
    
    # Now we have to "clip" these files.. (as same lenght R1 and R2 are error in bowtie..)
    echo "perl ${PerlHelpersPath}/trim1base3prime.pl extendedFrags.fastq" >&2
    perl ${PerlHelpersPath}/trim1base3prime.pl extendedFrags.fastq
    rm -f extendedFrags.fastq
    echo "perl ${PerlHelpersPath}/trim1base3prime.pl extendedFrags_RevCompl.fastq"  >&2
    perl ${PerlHelpersPath}/trim1base3prime.pl extendedFrags_RevCompl.fastq
    rm -f extendedFrags_RevCompl.fastq
    # returned files extendedFrags_3primeTrimmed.fastq extendedFrags_RevCompl_3primeTrimmed.fastq
    mv -f extendedFrags_3primeTrimmed.fastq FLASHED_1.fastq
    mv -f extendedFrags_RevCompl_3primeTrimmed.fastq FLASHED_2.fastq
    
    else
    
    printThis="Flash could not combine any paired end reads."
    printToLogFile
    
    #------------------------------------------
    echo '# end of "if [ -s extendedFrags.fastq ] ; then"'
    fi
    #------------------------------------------
    
    # We exit flashing with files :
    
    # Combined originally-mapping (and after-trimming-mapping if trimming was done) :
    # MAPPED_1.fastq MAPPED_2.fastq
    
    # Combined not-yet-mapping : (if trimming was done before flashing)
    # UNMAPPED_1.fastq UNMAPPED_2.fastq
    
    # Flashed reads :
    # FLASHED_1.fastq FLASHED_2.fastq
    
    echo
    
    #------------------------------------------
    echo '# end of "no unmapped reads" if then'
    fi

    else
    printThis="No FLASHing of unmapped paired end reads requested - skipping 'combining overlapping short reads' in flash !"
    printToLogFile  
    echo
    #------------------------------------------
    # end of if "No flashing"
    fi
    #------------------------------------------
    # end of if single end run
    fi
    #------------------------------------------
    
    # Now we have trimmed and/or flashed.
    # We have some of these files :
    
    # Now, we have at least SOME of these files :

    # MAPPED_1.fastq MAPPED_2.fastq (Combined originally-mapping and after-trimming-mapping if trimming was done)
    # M_FILTERED_1.fastq M_FILTERED_2.fastq (combined m-filtered reads)
    # UNMAPPED_1.fastq UNMAPPED_2.fastq (Combined not-mapping-not-flashable)
    # FLASHED_1.fastq FLASHED_2.fastq (Flash-extended reads)
    
    
    # Making catenation series (the files have to go in right order, so using * when combining is dangerous)
    
    mapped1=""
    mapped2=""
    unmapped1=""
    unmapped2=""
    filtered1=""
    filtered2=""
    flashed1=""
    flashed2=""
    
    # Making fastqc for the trimmed / flashed reads - to monitor what was done..
    if [ -s UNMAPPED_1.fastq ] ; then
    printThis="Running FastQC for unmapped reads.."
    printToLogFile
    printThis="${PipePath}/QC_and_Trimming.sh --fastqc --single ${SINGLE_END} --basenameR1 UNMAPPED_1 --basenameR2 UNMAPPED_2"
    printToLogFile
    ${PipePath}/QC_and_Trimming.sh --fastqc --single ${SINGLE_END} --basenameR1 UNMAPPED_1 --basenameR2 UNMAPPED_2
    
    if [ -s UNMAPPED_1_fastqc.html ] ; then
    mkdir UNMAPPED_1_fastqc
    mv -f UNMAPPED_1_fastqc.html UNMAPPED_1_fastqc/fastqc_report.html
    fi
    
    if [ -s UNMAPPED_2_fastqc.html ] ; then
    mkdir UNMAPPED_2_fastqc
    mv -f UNMAPPED_2_fastqc.html UNMAPPED_2_fastqc/fastqc_report.html
    fi
    
    fi
    
    if [ -s FLASHED_1.fastq ] ; then
    printThis="Running FastQC for flashed reads.."
    printToLogFile
    printThis="${PipePath}/QC_and_Trimming.sh --fastqc --single 1 --basenameR1 FLASHED_1"
    printToLogFile
    ${PipePath}/QC_and_Trimming.sh --fastqc --single 1 --basenameR1 FLASHED_1
    
    mkdir FLASHED_1_fastqc
    mv -f FLASHED_1_fastqc.html FLASHED_1_fastqc/fastqc_report.html  
    
    fi
   
    
    # Assuming that if we have R1, we have equally many in R2..
    if [ -s MAPPED_1.fastq ] ; then
        mapped1="MAPPED_1.fastq"
        mapped2="MAPPED_2.fastq"  
    fi
    
    if [ -s UNMAPPED_1.fastq ] ; then
        unmapped1="UNMAPPED_1.fastq"
        unmapped2="UNMAPPED_2.fastq"  
    fi
    
    if [ -s M_FILTERED_1.fastq ] ; then
        filtered1="M_FILTERED_1.fastq"
        filtered2="M_FILTERED_2.fastq"  
    fi

    if [ -s FLASHED_1.fastq ] ; then
        flashed1="FLASHED_1.fastq"
        flashed2="FLASHED_2.fastq"  
    fi

    
    cat ${mapped1} ${unmapped1} ${filtered1} ${flashed1} > READ1.fastq
    testedFile="READ1.fastq"
    doTempFileTesting
    rm -f *_1.fastq
    
    if [ "$SINGLE_END" -eq 0 ] ; then
    cat ${mapped2} ${unmapped2} ${filtered2} ${flashed2} > READ2.fastq
    testedFile="READ2.fastq"
    doTempFileTesting
    rm -f *_2.fastq
    fi
    
    cat READ1.fastq | grep -v "^$" > temp.fastq
    mv -f temp.fastq READ1.fastq
    
    if [ "$SINGLE_END" -eq 0 ] ; then
    cat READ2.fastq | grep -v "^$" > temp.fastq
    mv -f temp.fastq READ2.fastq
    fi

    echo
    echo "Read counts in R1 and R2 after trimming and/or flashing :"
    
    combo1count=0
    combo2count=0

    combo1count=$(( $( grep -c "" READ1.fastq )/4 ))
    if [ "$SINGLE_END" -eq 0 ] ; then
    combo2count=$(( $( grep -c "" READ2.fastq )/4 ))
    fi  
    
    echo "READ1.fastq : ${combo1count}"
    if [ "$SINGLE_END" -eq 0 ] ; then
    echo "READ2.fastq : ${combo2count}"
    fi

}


