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

setPathsForPipe(){
    
# #############################################################################

# This is the CONFIGURATION FILE to load in the needed toolkits ( conf/loadNeededTools.sh )

# #############################################################################

# Setting the needed programs to path.

setToolLocations=1
# If all the needed toolkits are already present in PATH, you can turn this feature off by setting :
# setToolLocations=0

# -----------------------------------------
# SETTING THE LOCATIONS OF THE TOOLKITS
# ------------------------------------------

if [ "${setToolLocations}" -eq 1 ]; then

# This can be done EITHER via module system, or via EXPORTING them to the path.
# If exporting to the path - the script does not check already existing conflicting programs (which may contain executable with same names as these)

# Change this to "0" if you want to use direct paths.
useModuleSystem=1
# useModuleSystem=1 : load via module system (set module names below)
# useModuleSystem=0 : load via direct paths (set paths below)

# ------------------------------------------
# PATHS_LOADED_AS_MODULES

# This is done, if you have set : useModuleSystem=1

if [ "${useModuleSystem}" -eq 1 ]; then

module purge
# Removing all already-loaded modules to start from clean table

module load samtools/0.1.19
# Supports all samtools versions in 0.* series. Does not support samtools/1.0 or above.

module load bowtie/1.0.0
# Supports all bowtie1 versions 1.* and 0.*
module load bowtie2/2.1.0
# Supports all bowtie2 versions

module load bedtools/2.17.0
# Supports bedtools versions 2.1* . Does not support bedtools versions 2.2*

module load ucsctools/1.0
# Supports ucsctools versions 1.* . Not known if supports also ucsctools versions 2.* (most probably not)
# The needed tools are :
#    bedGraphToBigWig
#    bedClip
#    bedGraphPack
#    bedToBigBed
#    bigBedToBed

module load flash/1.2.8
# Not known if would support other flash versions. Most probably will support.

module load fastqc/0.11.4
# Will not support fastqc versions 0.10.* or older

module load trim_galore/0.3.1
# Not known if would support other trim_galore versions. Most probably will support.


# module load cutadapt/1.2.1
# If your trim_galore module does not automatically load the needed cutadapt,
# uncomment this line

# Not known if would support other cutadapt versions. Most probably will support.

module load perl/5.18.1
# Most probably will run with any perl

module list
perl --version

else

# ------------------------------------------
# EXPORT_PATHS_WITHOUT_MODULE_SYSTEM 

# This is done, if you have set : useModuleSystem=0

# Note !!!!!
# - the script does not check already existing conflicting programs within $PATH (which may contain executable with same names as these)

export PATH=$PATH:/usr/local/pkgbin/samtools/0.1.19/bin
export PATH=$PATH:/usr/local/pkgbin/bowtie/1.0.0/bin
export PATH=$PATH:/usr/local/pkgbin/bowtie2/2.1.0/bin
export PATH=$PATH:/usr/local/pkgbin/bedtools/2.17.0/bin
export PATH=$PATH:/usr/local/pkgbin/ucsctools/1.0/bin
export PATH=$PATH:/usr/local/pkgbin/flash/1.2.8/bin
export PATH=$PATH:/usr/local/pkgbin/fastqc/0.11.4/bin
export PATH=$PATH:/usr/local/pkgbin/trim_galore/0.3.1/bin
export PATH=$PATH:/usr/local/pkgbin/cutadapt/1.2.1/bin
export PATH=$PATH:/usr/local/pkgbin/perl/5.18.1/bin

# See notes of SUPPORTED VERSIONS above !

echo $PATH
perl --version

fi
fi

}

