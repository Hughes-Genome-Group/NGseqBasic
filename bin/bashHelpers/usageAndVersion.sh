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


version(){
    nameOfScript="VERSION 20 FirstPortable"
    manualFileName="DNasePipeUserManual_VS_101_090415.pdf"
    manualURLpath="http://userweb.molbiol.ox.ac.uk/public/telenius/MANUAL_for_pipe_030214/DNasePipeUserManual_VS_101_090415.pdf"
    pipesiteAddress="http://userweb.molbiol.ox.ac.uk/public/telenius/PipeSite.html"
    
    versionInfo="\n${nameOfScript} NGseqBasic.sh\nmanual available in ${manualURLpath}\nUpdates (bug fixes etc) listed in : ${pipesiteAddress}\n"
    echo -e "${versionInfo}" > versionInfo.txt
    echo "<p>${nameOfScript} NGseqBasic.sh </p>" > versionInfoHTML.txt
    echo "<p> User manual - to understand the pipeline and the output :  <a href=\"${manualURLpath}\" >${manualFileName}</a> </p>" >> versionInfoHTML.txt
    echo "<p> Updates (bug fixes etc) listed in : <a href=\"${pipesiteAddress}\" >${pipesiteAddress}</a></p>" >> versionInfoHTML.txt
}

usage(){
    
    version
    echo -e ${versionInfo}
    rm -f versionInfoHTML.txt
    rm -f versionInfo.txt

echo "FOR ATAC/DNASE/CHIP SEQUENCING DATA"
echo
echo "[1](bam-to-fastq) --> [1b] fastqc --> (trimming, fastqc again) --> bowtie --> bowtie_out.bam, filtered.bam (no blacklisted*, no dupl, proper pairs) --> bedgraph, bigwig --> [2](data hub)"
echo "Data can enter as bam [1], and/or as fastq [1b]. Step [2] is optional."
echo "Note that the generated bigwig and bedgraph files contain SINGLE INTERVAL for each properly paired READ PAIR (mimicking the original DNaseI cut fragments)"
echo "*) Blacklisted region filtering works only for "
for g in $( seq 0 $((${#genomesWhichHaveBlacklist[@]}-1)) ); do echo -n "${genomesWhichHaveBlacklist[$g]} "; done
echo "at the moment"
echo
echo "Is to be ran in the command line like this :"
echo "qsub -cwd -o qsub.out -e qsub.err -N MyPipeRun < ./run.sh"
echo "Where run.sh is oneliner like : '${PipePath}/DnaseAndChip_pipe_1.sh --genomes mm9,mm10' "
echo -n "Supported genomes : "
for g in $( seq 0 $((${#supportedGenomes[@]}-1)) ); do echo -n "${supportedGenomes[$g]} "; done
echo "(see * note above about blacklisted region filtering)"
echo
echo "Run the script in an empty folder - it will generate all the files and folders it needs."
echo
echo "THE SCRIPT NEEDS 1-3 PARAMETER FILES :"
echo "Have a look at the Example set of parameter files with this command :"
echo "ls /home/molhaem2/telenius/Jelena_DNase_pipe/ExampleRunSetup/*"
echo "OBLIGATORY : either PIPE_bamPaths.txt or PIPE_fastqPaths.txt (or both, see below), or PIPE_mappedBamPaths.txt if you start with already bowtie-mapped bam files"
echo "OPTIONAL : PIPE_hubbing.txt, PIPE_hubbingSymbolic.txt"
echo "RERUN : if you want to rerun the Footprint or Peak call or Hubbing, you need to provide PIPE_previousRunPaths.txt (as the OBLIGATORY file) - start also this run in empty folder (not in the folder of the previous run)"
echo
echo "Examples to create parameter files - HANDS ON TUTORIAL : http://userweb.molbiol.ox.ac.uk/public/telenius/MANUAL_for_pipe_030214/DnaseCHIPpipe_TUTORIAL.pdf"
echo
echo "OPTIONAL FLAGS FOR TUNING THE PIPE RUN :"
echo
#echo "DnaseAndChip_pipe_1.sh --genomes mm9,mm10 -h --help --bowtie1 --bowtie2 --noBowtie --onlyHub --chip --dnase --noUnpair --unpair --trim --noTrim --maxins 100 --mergeFP 1 --contigFP 2 --depthFP 3 --window 3 --incr 1"
#echo "DnaseAndChip_pipe_1.sh --genomes mm9 -h --help --noBowtie --onlyHub --noUnpair --unpair --trim --noTrim --sam2bw --maxins 300 --mergeFP 100 --contigFP 10 --depthFP 10 --window 9 --incr 3 --mergeP 10 --contigP 70 --depthP 20"
echo "DnaseAndChip_pipe_1.sh --genomes mm9 -h --help --chunkmb 256 -p 3 --lanes 1 --noBowtie --onlyPeakCall --onlyFPandPC --saveUnpaired --saveUnfiltered --trim --noTrim --flash --noFlash --noSam2bw --maxins 300 --mergeFP 100 --contigFP 10 --depthFP 10 --window 9 --incr 3 --mergeP 0 --contigP 20 --depthP 40 --normDepth"
echo
echo "Default behavior :"
echo
echo "DnaseAndChip_pipe_1.sh --lanes 1 --chunkmb 256 -m 2 --trim --flash --blacklistFilter --windowSize 300 --windowIncr 30 --maxins 350 --mergeFP 100 --contigFP 10 --depthFP 10 --windowFP 1 --mergeP 0 --contigP 20 --depthP 40"
echo
echo "HELP"
echo "-h, --help : prints this help"
echo
echo "RUN COMMANDS"
echo "Example queue run :"
echo "   qsub run.sh -o qsub.out -e qsub.err , where run.sh is oneliner 'NGseqBasic.sh --genomes mm9' "
echo "Example nohup run :"
echo "   nohup NGseqBasic.sh --genomes mm9 1> qsub.out 2> qsub.err & "
echo "More detailed run commands in hands-on-tutorial : http://userweb.molbiol.ox.ac.uk/public/telenius/MANUAL_for_pipe_030214/DnaseCHIPpipe_TUTORIAL.pdf "
echo
echo "OUTPUT LOG FILE NAMES"
echo "--outfile qsub.out (the STDOUT log file name in your RUN COMMAND - see above )"
echo "--errfile qsub.err (the STDERR log file name in your RUN COMMAND - see above )"
echo ""
echo "ACCESSIBILITY SETTINGS"
echo "--orangeBlue (use orange-blue color scheme instead of red-green)"
echo "--redGreen   (use the OLD redGreen colors - very colorblind unfriendly, i.e. the default colors of runs before 2019)"
echo "Default : not use either of above (use pink-green colors, which are colorblind friendly)"
echo ""
echo "FASTQ SETTINGS"
echo "--gz (input files are provided in file.fastq.gz compressed format )"
echo "--lanes 2 (set this to be the number of lanes, if there are more than 1 lanes in your fastq files)"
echo "--singleEnd - to run single end sequencing files (default behavior is paired end files)"
echo ""
echo "BOWTIE SETTINGS"
echo "-p 3 : to how many processors we will parallelise the bowtie part of the run"
echo "--bowtie1 / --bowtie2 (default is bowtie1 - decide if bowtie1 or bowtie2 is to be used. bowtie2 is better to long reads - read lenght more than 70b, fragment lenght more than 350b)"
echo "--singleEnd - to run single end sequencing files (default behavior is paired end files)"
echo "--chunkmb "${BOWTIEMEMORY}" - memory allocated to Bowtie, defaults to 256mb - only affects bowtie1 run"
echo "-M 2 run with bowtie parameter M=2 (if maps more than M times, report one alignment in random) - only affects bowtie1 run"
echo "-m 2 run with bowtie parameter m=2 (if maps more than m times, do not report any alignments) - only affects bowtie1 run"
echo "-m and -M are mutually exclusive."
echo "--maxins 350 : sets the TRUE fragment (max) lenght to 350bp. Bowtie1 default 250, bowtie2 default 500, this script default 350. (Bowtie mappings resulting in fragments longer than --maxins are not reported)"
echo "-v 3 : allow up-to-this-many total mismatches per read (ignore base qualities for these mismatches). "
echo "       cannot be combined to --seedlen, --seedmms or --maqerr (below)."
echo "--seedlen 28 - alignment seed lenght (minimum 5 bases) . Seed is the high-quality bases in the 5' end of the read. Default 28 (bowtie1), 20 (bowtie2)."
echo "--seedmms 2 - max mismatches within the seed (see the 'seed' above). Allowed 0,1,2,3 mismatches in bowtie1 - default 2, allowed 0,1 in bowtie2 (per each multi-seed alignment) - default 0. "
echo "--maqerr 70 - only in Bowtie1 - max total quality values at all mismatched read positions throughout the entire alignment (not just in seed)"
echo ""
echo "TURN FEATURES ON/OFF"
echo "--noBowtie runs only pipe AFTER bowtie (assumes PIPE_mappedBamPaths.txt, see above)"
echo "--flash**/noFlash (run/do-not-run Flash - for unmapped reads 'try to merge overlapping short reads to longer single end read' to enhance mapping)"
echo "**) NOTE : combination --flash --noTrim is not recommended (can result in combining reads on the sites of ADAPTERS instead of the reads themselves). Use with caution."
echo -n "--blacklistFilter/noBlacklistFilter - supported genomes "
for g in $( seq 0 $((${#genomesWhichHaveBlacklist[@]}-1)) ); do echo -n "${genomesWhichHaveBlacklist[$g]} "; done
echo "(filter/do-not-filter blacklisted regions out from the final bam file)"
echo "--footPrint : generate also footprint tracks (see footprint(FP) parameters below)"
echo "--peakCall : generate also peak call (see peak call(P) parameters below)"
echo ""
echo "ADAPTER TRIMMING SETTINGS (for reads which don't map without trimming)"
echo "--nextera  : use Nextera adaptors in trimming, instead of standard illumina PE adapters"
echo "--trim/noTrim** (run/do-not-run TrimGalore for the data - Illumina PE standard adapter filter, trims on 3' end)"
echo "**) NOTE : combination --flash --noTrim is not recommended (can result in combining reads on the sites of ADAPTERS instead of the reads themselves). Use with caution."
echo "--trim5 (run trimming also for 5' end of the data - Illumina PE standard adapter filter, combination of TrimGalore and cutadapt) "
echo "--ada3read1 SEQUENCE --ada3read2 SEQUENCE  : custom adapters 3' trimming, R1 and R2 (give both if PE custom trimming is needed, SE trimming needs only R1) - these adapters will be used instead of Illumina default / atac adapters. SEQUENCE has to be in CAPITAL letters ATCG"
echo "--ada5read1 SEQUENCE --ada5read2 SEQUENCE  : custom adapters 5' trimming, R1 and R2 (give both if PE custom trimming is needed, SE trimming needs only R1) - these adapters will be used instead of Illumina default / atac adapters. SEQUENCE has to be in CAPITAL letters ATCG"
echo ""
echo "WINDOWING IN TRACK VISUALISATION"
echo "Default windowing is 300b window and 10% increment (see --windowSize and --windowIncr parameters to set non-default window)"
echo "--noWindow : no windowing - instead plot 1b resolution raw read depths"
echo "--windowSize 300 : custom window size (instead of default 300b) - this value has to be even number (or is rounded into one)."
echo "                  If no custom INCREMENT is given, increment is set to 10% of windowSize , i.e. the windowSize value has to be divisible by 10 (or will be rounded to behave as such)."
echo "--windowIncr 30 : custom window increment (instead of default 10%). The value of windowIncr should be even number, given in BASES (10 bases : windowIncr=10), not percentages."
echo "                  If no custom windowSize is given, the windowSize will be set to 300bases, and increment rounded so that 300bases is divisible by the given increment (allows only values 10 and 30)."
echo "Note that custom --windowSize and --windowIncr should be compatible - both need to be even numbers, and windowSize should be exact multiply of windowIncr (if not - will be rounded to reach those requirements)"
echo ""
echo "SAVE EXTRA BAM / FASTQ / BDG,WIG FILES - default output is : duplicate(+blacklisted) -filtered bam, including ONLY mapped reads (single-end), or PROPER PAIRS (paired end)"
echo "--saveUnmapped (save fastq files UNMAPPED_1.fastq UNMAPPED_2.fastq for not-paired-mapped reads : only for PE data)"
echo "--saveUnpaired (save bam files singleEnd_bowtie_READ1.bam singleEnd_bowtie_READ2.bam for unpaired mapped reads - not blacklisted/duplicate-filtered)"
echo "--saveUnpairedFiltered (save bam files singleEnd_bowtie_filtered_READ1.bam singleEnd_bowtie_filtered_READ1.bam for unpaired mapped reads - blacklisted/duplicate-filtered)"
echo "--saveUnfiltered (save original bam file bowtie_out.bam - straight from bowtie output)"
echo "--saveUnfilteredMapped (save original bam file mapped.bam - straight from bowtie output. Do not print out unmapped reads (single end) / non-proper pairs (paired end).)"
echo "--saveUntrimmed (run EXTRA bowtie-run before trimming, save the bam file UNTRIMMED_bowtie_out.bam - straight from before-trimming bowtie output - includes all reads)"
echo "--saveUntrimmedMapped (run EXTRA bowtie-run before trimming, save the bam file UNTRIMMED_onlyMapped_bowtie_out.bam - straight from before-trimming bowtie output. Do not print out unmapped reads (single end) / non-proper pairs (paired end).)"
echo "--saveBDG (save unpacked BDG files - default : save only BIGWIG packed versions of these files"
echo ""
echo "RERUN SETTINGS"
echo
echo "--pyramidRerun : if you are running peak calls / footprints via PYRAMID pipeline, set this on."
echo "--onlyPeakCall : only run Peak Call (to finetune peak call part) - expects PIPE_previousRunPaths.txt (see above)"
echo "--onlyFPandPC : only run Peak Call and FootPrinting (to finetune FootPrint windowing and peak call parameters) - expects PIPE_previousRunPaths.txt (see above)"
echo "--onlyHub : only generates data hub - NOTE !! this DOES NOT include data from any 'footprint / peak call rerun' folders - it hubs only the ORIGINAL data set - (expects PIPE_previousRunPaths.txt, see above)"
echo 
echo "FOOTPRINT parameters (window,increment) :"
echo
echo "--windowFP 1  : (windowing the signal for visualisation) - sliding window, where overlap 2*WINDOW, and window size 4*WINDOW, resulting in 2*WINDOW size increments in graph. Value 0 means NO WINDOWING (only raw 1b resolution track is given)"
echo
echo "PEAK CALL parameters (mergeP,contigP,depthP) :"
echo
echo "--depthP 40  : (first filter)  RANGE of depthP - Each region should have at least 40 reads to be included to the peak call. "
#echo "--normDepth  : (first filter)  If --normDepth is stated, it will interpret depthP = 'reads / million reads in sample'."
echo "--mergeP 0   : (second filter) Merges regions separated by 0bp or less."
echo "--contigP 20 : (third filter)  Regions (merged, depth-filtered) should be at least 20bp wide to be considered."
echo
echo "More info : hands-on tutorial : http://userweb.molbiol.ox.ac.uk/public/telenius/MANUAL_for_pipe_030214/DnaseCHIPpipe_TUTORIAL.pdf, comprehensive user manual : http://userweb.molbiol.ox.ac.uk/public/telenius/MANUAL_for_pipe_030214/DNasePipeUserManual_VS_100_180215.pdf , and comment lines (after the subroutine descriptions) in the script ${PipePath}/DnaseAndChip_pipe_1.sh"
echo

exit 0

}

parameterUsage(){

    version
    echo -e versionInfo
    rm -f versionInfoHTML.txt

echo "OBLIGATORY : either PIPE_bamPaths.txt or PIPE_fastqPaths.txt (or both, see below), or PIPE_mappedBamPaths.txt if you start with already bowtie-mapped bam files"
echo "OPTIONAL : PIPE_hubbing.txt, PIPE_hubbingSymbolic.txt"
echo "RERUN : if you want to rerun the Footprint or Peak call or Hubbing, you need to provide PIPE_previousRunPaths.txt"
echo
echo "----------------------------------------------------------------------------------"
echo "PIPE_bamPaths.txt"
echo "SampleName <TAB> /file/path/for/sample.bam -One line per sample, end with newline"
echo "----------------------------------------------------------------------------------"
echo "PIPE_fastqPaths.txt"
echo "SampleName <TAB> /fastq/file/path/for/read1.fastq  <TAB> /fastq/file/path/for/read2.fastq -One line per sample, end with newline"
echo "If you have more than 1 lane in sequencing (=more than 2 fastq files), use --lanes flag (below) and write PIPE_fastqPaths.txt like this (for 2 lanes):"
echo "SampleName <TAB> /path/read1_lane1.fastq,/path/read1_lane2.fastq <TAB> /path/read2_lane1.fastq,/path/read2_lane2.fastq"
echo "----------------------------------------------------------------------------------"
echo " --> provide either PIPE_bamPaths.txt or PIPE_fastqPaths.txt, or both parameter files (if some your samples have bams and others fastqs as starting files)"
echo "----------------------------------------------------------------------------------"
echo "PIPE_mappedBamPaths.txt"
echo "SampleName <TAB> /file/path/for/sample.bam <TAB> genomeBuild -One line per sample, end with newline (genomeBuild = build to which the data is mapped to)"
echo "----------------------------------------------------------------------------------"
echo "PIPE_previousRunPaths.txt"
echo "SampleName <TAB> /folder/path/of/previous/run genomeBuild -One line per sample, end with newline (genomeBuild = build to which the data is mapped to)"
echo "This file path is the path to the run FOLDER generated in previous run : /my/previous/running/dir/SAMPLE_DIR "
echo "----------------------------------------------------------------------------------"
echo "If UCSC data hub is to be created - provide either PIPE_hubbing.txt or PIPE_hubbingSymbolic.txt (see below)"
echo "In order to use PIPE_hubbingSymbolic.txt - ensure you can provide an EVER-LASTING file path - the symbolic links will break if you move the data away from this path"
echo "The UCSC data hub can be generated also afterwards by running this same script with flag --onlyHub (see above)"
echo "----------------------------------------------------------------------------------"
echo "PIPE_hubbing.txt - optional oneliner (end with newline), if UCSC data hub is to be created."
echo "hubName <TAB> /public/file/path/for/data/hub/folder <TAB> /public/file/path/to/store/visualisation/bigwig/files"
echo "As it is advisable to have the hub and bigwigs in separate folders, these 2 cannot point to the same folder"
echo "(Permanent location for bigwigs = not so easy to break the data hub. Freedom to move around and reuse the hub = more flexible hub.)"
echo "----------------------------------------------------------------------------------"
echo "PIPE_hubbingSymbolic.txt - optional oneliner (end with newline), if UCSC data hub is to be created"
echo "hubName <TAB> /public/file/path/for/data/hub/folder <TAB> /public/file/path/to/store/visualisation/files/and/symbolic/links <TAB> /file/path/to/PERMANENT/BIGWIG/STORAGE/FOLDER/yourHubName"
echo "The permanent storage folder can be located anywhere in the CBRG system - as long as it is visible to the /public/ directory (any /hts/dataX/username area dues fine). "
echo "Do not MOVE the data stored in your PERMANENT/STORAGE/FOLDER - your symbolic links will stop functioning and break the data hub - if you are moved from (say) data1 to data3 - leave this folder behind !"
echo "----------------------------------------------------------------------------------"
echo
echo "HANDS ON TUTORIAL : http://userweb.molbiol.ox.ac.uk/public/telenius/MANUAL_for_pipe_030214/DnaseCHIPpipe_TUTORIAL.pdf"
echo

exit 0

}


