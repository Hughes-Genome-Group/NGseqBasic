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

=head1 NAME

   fastq_scores.pl 

=head1 SYNOPSIS

   fastq_scores.pl -i <input file> -r <num reads to examine>

=head1 EXAMPLE

   fastq_scores.pl -i mouse.txt -r 25000

=head1 OPTIONS

   -i input file: sequence file (FASTQ format)
   -r num reads to examine (default: unlimited)
 
=head1 DESCRIPTION

   This script will attempt to determine the scoring scheme used in a given fastq file

=head1 AUTHOR

   (c) S. McGowan, Computational Biology Research Group, Oxford
   www.cbrg.ox.ac.uk

=head1 Update Record

   Aug 2010  001  S.McGowan  first written
   May 2011  002  S.McGowan  added option to stop after defined num reads
   Sep 2011  003  S.McGowan  added Phred+64 scores (but not tested on genuine Phred+64 file yet!)
   Feb 2014       J Telenius -modified output to fit to the DnaseAndChip_pipe1.sh pipeline

=cut

#-------------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;

$| = 1;

my $help=0;
my $man=0;

GetOptions(
"i=s"=>\my $seq_reads_file,	
"r:i"=>\my $read_limit,	
'h|help'=>\$help,
'man'=>\$man
); 

# Printing version name
# print "";
# print "Perl version :";
# print "$^V";
# print "";

pod2usage(1) if $help;
pod2usage(-verbose=>2) if $man;
pod2usage(2) unless ($seq_reads_file);


unless(-e $seq_reads_file)
{
#	print "\nUnknown file: $seq_reads_file\nUSAGE: perl /package/cbrg/lib/fastq_scores.pl /full/path/to/fastq_file.fq\n\n";
	exit();
}

#print "\n\tThis script will attempt to determine the\n\tscoring scheme used in a given fastq file:\n";
#print "\n";

my $line_count = 1;
my $total_count = 0;

my $min_score = 999;
my $max_score = 0;

open (FASTQ, "$seq_reads_file") or die "Couldn't open file $seq_reads_file\n\n"; 
while (<FASTQ>)
{
	my $line = $_; 	
	chomp $line;			
	if ($line_count < 4)
	{
		$line_count ++;
	}	
	else        
	{
		$total_count ++;
		my @indiv_scores = split(//,$line);
		foreach my $ascii_score (@indiv_scores)
		{
			my $score_value = &decode_scores($ascii_score);
			if ($score_value < $min_score){$min_score = $score_value;}
			if ($score_value > $max_score){$max_score = $score_value;}
		}
		$line_count = 1;
		
		if ($read_limit)
		{
			if ($total_count >= $read_limit)
			{
				my $prob_scheme = &get_scheme($min_score, $max_score, $total_count);
				#print "\tReads processed:$total_count | Min score:" . chr($min_score) . " ($min_score)  | Max score:" . chr($max_score) . " ($max_score) $prob_scheme\r";
				print "$prob_scheme";
				last;
			}
		}
		elsif ( ($total_count > 25000)  and (($total_count % 1000) == 0) )
		{
			#my $prob_scheme = &get_scheme($min_score, $max_score, $total_count);
			#print "\tReads processed:$total_count | Min score:" . chr($min_score) . " ($min_score)  | Max score:" . chr($max_score) . " ($max_score) $prob_scheme\r";
			#print "$prob_scheme";
		}
	}
}
close FASTQ;

#print "\n\n";

exit;
#---------------------------------


sub decode_scores
{
	my ($ascii_score) = @_;
	return (ord($ascii_score));
}

sub get_scheme
{
		my ($min_score, $max_score, $total_count) = @_;
		
		#if ($max_score <= 73) {return "| sanger (Phred+33)."}
		#elsif ($max_score <= 74) {return "| Illumina 1.8+ (Phred+33)."}
		#elsif (($min_score >= 58) and ($min_score < 64)) {return "| solexa."}
		#elsif (($min_score == 64) or ($min_score == 65)) {return "| Illumina (Phred+64) (v1.3+)"}
		#elsif ($min_score >= 66)  {return "| Illumina (Phred+64)"}
		#else {return "| who knows what scheme is being used here!"}
		
		if ($max_score <= 73) {return "--phred33-quals"}
		elsif ($max_score <= 74) {return "--phred33-quals"}
		elsif (($min_score >= 58) and ($min_score < 64)) {return "--solexa-quals"}
		elsif (($min_score == 64) or ($min_score == 65)) {return "--phred64-quals"}
		elsif ($min_score >= 66)  {return "--phred64-quals"}
		else {return "--phred33-quals"}
}



