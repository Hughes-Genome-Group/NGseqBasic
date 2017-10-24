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

use Pod::Usage;
use Getopt::Long;
use strict;

my $window_size;
my $window_incr;
my %chrlength;
my $username;
my $track_name;
my $chr_lengths_file;

# This script is based on sam2bwPE.pl by Jim Hughes (using version which was available in CBRG system in /package/cbrg/bin/sam2bwPE.pl at 17Apr2015 )
# Modifications - to make it read in bed file instead of sam, and only do windowing and no other fancy stuff, by Jelena Telenius, April 2015

&GetOptions (
	     "bed=s"=>\my $bedfile,
	     "window=i"=> \$window_size,
	     "name=s"=>\$track_name,
	     "genome=s"=>\$chr_lengths_file,
	     "inc=i" => \$window_incr);

# Printing version name
print "";
print "Perl version :";
print "$^V";
print "";
    
# set some defaults:
unless($window_size) {$window_size = 300;}
unless($window_incr) {$window_incr = 30;}
unless($track_name) {$track_name = 'track_name';}

open (SIZES, $chr_lengths_file);
	
while (<SIZES>)
{
	chomp;
	my ($gchr, $gsize)= split(/\t+/);
	print "$gchr $gsize\n";

	$gchr =~ s/chr//gi;
	unless ($gchr =~ /M/gi)  { $chrlength{$gchr}=$gsize; }		
}	
close SIZES;
	
my $pos_fix = int($window_size/2);

open (OUTPUT, ">$track_name.windowed.bed");

open(INFO1, $bedfile) || die $!;

my %bins;
my %read_posns;
my $bintotal;

while (<INFO1>) {
	&HASHIT;
}

close INFO1;

my $binsnumbers;

foreach my $dChr(sort by_number keys %bins) {
	foreach my $dcoor(sort by_number keys %{$bins{$dChr}}) {
		
		my $dcount = $bins{$dChr}{$dcoor};
		my $chromolength = $chrlength{$dChr};

		my $adjustcoor = $dcoor + $pos_fix;   ####move coordinate from start of window to middle of window
		my $adjustcoorend = $adjustcoor + $window_incr;

		unless ($adjustcoor > $chromolength){
			print OUTPUT "chr$dChr\t$adjustcoor\t$adjustcoorend\t$dcount\n";
			$binsnumbers++;
			$bintotal += $dcount;
			}	
}}

my $avwindowSig = 0;
if ($binsnumbers) {
     my $avwindowSig = $bintotal/$binsnumbers;
}

print "Average window signal is $avwindowSig\n";


##########subs#########


sub HASHIT {

	my ($chrbed, $startbed, $endbed)=split(/\s+/);
	if ($chrbed =~ /_/g){next;}
	$chrbed =~ s/chr//;
	
	unless(exists $chrlength{$chrbed}){next;}

	
	my $int = int($startbed / $window_size);
	my $start_bin = ($int * $window_size);
	my $diff = $startbed - $start_bin;
	my $incr = (int(($window_size - $diff) / $window_incr) * $window_incr);
	$start_bin -= $incr;
	
	for (my $bin=$start_bin; $bin<($endbed); $bin+=$window_incr){
	#for (my $bin=$start_bin; $bin<($start_bin+$window_size); $bin+=$window_incr){
		#unless (($chrbed =~ /M|m/)||($bin <= 0)) {
		unless ($bin <= 0) {
		$bins{$chrbed}{$bin} ++;
		}
	}

}

sub by_number {
	($a <=> $b);
	}
