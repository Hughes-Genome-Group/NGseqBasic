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

use strict;

# Jelena Telenius 2014 - Clipping the last base from 3' end
# James Davies 2014 - Reading in the fastq in batches of four lines (from script dpnII2E.pl )

# Printing version name
print "";
print "Perl version :";
print "$^V";
print "";

my $line_counter = 0; my $f_counter = 0; my %hash; my $flag;
my @line_labels = qw(name seq spare qscore);

my $filename = $ARGV[0];
$filename =~ /(.*)\.(fastq|fq)/; my $filename_out= $1."_3primeTrimmed.$2";

unless (open(FH, $filename)) {print "Cannot open file $filename\n"; exit;}

# opens a file in append modes for the output of the data
open FHOUT, ">$filename_out" or die $!;   
  
while ($hash{$line_labels[$f_counter]}=<FH>)  #assigns each fq line to the hash in batches of 4
{
chomp $hash{$line_labels[$f_counter]};
$f_counter++; $line_counter++;

if ($f_counter==4)
    {
   
         print FHOUT $hash{"name"}."\n";
         print FHOUT substr($hash{"seq"}, 0, -1)."\n";  
         print FHOUT "+\n";
         print FHOUT substr($hash{"qscore"}, 0, -1)."\n";  


    $f_counter=0

    }

}
