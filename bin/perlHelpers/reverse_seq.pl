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

# Printing version name
print "";
print "Perl version :";
print "$^V";
print "";

# Jelena Telenius 2014 - Reverse complementing of the reads (from Jelena's script sam2fastq_overwrite.pl - 2013)
# James Davies 2014 - Reading in the fastq in batches of four lines (from script dpnII2E.pl )

my $line_counter = 0; my $f_counter = 0; my %hash; my $flag;
my @line_labels = qw(name seq spare qscore);

my $filename = $ARGV[0];
$filename =~ /(.*)\.(fastq|fq)/; my $filename_out= $1."_RevCompl.$2";

unless (open(FH, $filename)) {print "Cannot open file $filename\n"; exit;}

# opens a file in append modes for the output of the data
open FHOUT, ">$filename_out" or die $!;   
  
while ($hash{$line_labels[$f_counter]}=<FH>)  #assigns each fq line to the hash in batches of 4
{
chomp $hash{$line_labels[$f_counter]};
$f_counter++; $line_counter++;

if ($f_counter==4)
    {

         # name formats @HISEQ2000:376:C2399ACXX:8:1101:1749:1893 1:N:0:GAGTTAGT run1   <-- we have only these (R1 reads output from FLASH)
         # name formats @HISEQ2000:376:C2399ACXX:8:1101:1749:1893 2:N:0:GAGTTAGT run2   <-- we want to output these (R2 read output)
         #           /(.*):(.*):(.*) :(.*):(.*):(.*):(\d++) (\d):(.*):(.*):(.*)/
         # Illumina name format from website : http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm
         
         # Here the names "1" - i.e. READ1 identifier, is changed to "2" - READ2 names.
         if ($hash{"name"} =~ /(@.*) (1):(.*)/)
             {
                  my $nameBegin = $1;
                  my $nameEnd = $3;
                  
                  $hash{"name"} = "".$nameBegin." 2:".$nameEnd;
             }
          
          
          
            #Reversing
	    my $seqRev = scalar reverse ($hash{"seq"}) ;
	    my $score = scalar reverse ($hash{"qscore"}) ;
	    
	    #Complementing
	    
	    #For each of the four bases - first changing the base to lower case complement, then changing back to upper case after complementing.
	    
	    #$seq =~ s/A/t/g ;
	    #$seq =~ s/T/a/g ;
	    #$seq =~ s/C/g/g ;
	    #$seq =~ s/G/c/g ;
	    
	    #$seq =~ s/t/T/g ;
	    #$seq =~ s/a/A/g ;
	    #$seq =~ s/g/G/g ;
	    #$seq =~ s/c/G/g ;
	    
	    #Faster version of previous lines :
	    
	    #$seq =~ s/A/t/g ;
	    #$seq =~ s/T/A/g ;
	    #$seq =~ s/t/T/g ;
	    
	    #$seq =~ s/C/g/g ;
	    #$seq =~ s/G/C/g ;
	    #$seq =~ s/g/G/g ;
	    
	    #Seems to be better to use tr which is "translitteration"
	    
	    $seqRev =~ tr/ATCG/atcg/ ; #to be on the safe side, making this in 2 steps
	    $seqRev =~ tr/atcg/TAGC/ ;   
    
         print FHOUT $hash{"name"}."\n";
         print FHOUT $seqRev."\n";  
         print FHOUT "+\n";
         print FHOUT $score."\n";  


    $f_counter=0

    }

}
