##########################################################################
# Copyright 2016, Jelena Telenius (jelena.telenius@imm.ox.ac.uk)         #
#                                                                        #
# This file is the README of NGseqBasic .                                #
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

Installation instructions :

1) Ask Jelena (jelena __ telenius __ at __ gmail __ com) to send you the tar.gz file including the codes.

2) Unpack with this command :
    tar --preserve-permissions -xzf NGseqBasic.tar.gz

3) Fill in the locations (or modules) of the needed tools (bowtie, fastqc etc) to the conf/config.sh file
    nano NGseqBasic/conf/config.sh       # Instructions as comment lines in the config.sh file

4) Fill in your server address to the conf/config.sh file
    nano NGseqBasic/conf/config.sh       # Instructions as comment lines in the config.sh file

5) Add the main script NGseqBasic.sh to your path or BASH profile (optional), f.ex :
    export PATH:${PATH}:/where/you/unpacked/it/NGseqBasic/NGseqBasic.sh

6) Start using the pipe ! (no installation needed)

7) Good place to start is the pipeline's help :
    NGseqBasic.sh --help

8) Below web site provides a test data set, hands-on tutorials, full manual, and other documentation !
   http://userweb.molbiol.ox.ac.uk/public/telenius/NGseqBasicManual/outHouseUsers/
   
9) Direct link to the test data set :
   http://userweb.molbiol.ox.ac.uk/public/telenius/NGseqBasicManual/outHouseUsers/testDataGzip_010117.tar.gz
   
   
