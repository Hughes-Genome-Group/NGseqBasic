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

use strict ;

# Printing version name
print "";
print "Perl version :";
print "$^V";
print "";

require Time::localtime ;
my $timestamp ;
my $output1 ;
my $output2 ;
my $output3 ;

my $counter = 0 ;
my $datalinecounter = 0 ;
my $printcount = 0 ;
my $orphancountR1 = 0 ;
my $orphancountR2 = 0 ;
my $duplicatecountR1 = 0 ;
my $duplicatecountR2 = 0 ;

my $waitingForR2 = 0 ;
my $waitingForR1 = 0 ;

my $dataR1inWaiting = "" ;
my $dataR2inWaiting = "" ;

my $name ;
my $data ;
my $flag ;

my $prevname ="";
my $pairprinted =0;
my $prevdata ;
my @prevline ;

my $READ1flag = hex("0x0040") ; # The bit which needs to be "1" if the read data line belongs to READ1
my $READ2flag = hex("0x0080") ; # The bit which needs to be "1" if the read data line belongs to READ2
my $REV_COMPLflag = hex("0x0010") ; # If this bit is "1", the data needs to be reverse complemented before printing

# This script assumes input to be name-wise SORTED sam file (or raw file from sequencing) - in any case so, that reads from same pair are grouped together.
# The reads do not need to be in order R1 R2 R1 R2, as long as same read names are grouped together.

{

$output1 = ">READ1.fastq";
$output2 = ">READ2.fastq";

}

open PRINTER1, $output1 ;
open PRINTER2, $output2 ;

print STDOUT "Beginning to read data..\n" ;

LINE : while (<>) {
#Saving data from input file, printing out the ones which are already in order.
    chomp;
    
    $counter++ ;
    if ( $counter % 1000000 == 0){
    print STDOUT "One million lines processed\n";    
    }
    
    my @thisline=split/\t/ ;
    
    unless ($thisline[0] =~ /^@/) {
        
        $datalinecounter++ ;
        
	# save old name
	$prevname = $name ;
	# set new name
        $name = "$thisline[0]" ;
        $flag = $thisline[1] ;
	
	#-----------------------------------------------------------------------------------------
	# When read name changes, we zero parameters.
	#-----------------------------------------------------------------------------------------
	
	# If this read was already printed, we will just move to next line (i.e. there was already a read pair red in with this name)
	# Here we of course assume input is sorted by name.
	
	# When read name changes.. - Zeroing counters and setting new readname.
	if ($prevname ne $name){
	    $pairprinted = 0;
	    
	    # Notifying orphan reads before zeroing the counter :
	    
	    if ($waitingForR1==1) {
		print STDOUT "Discarding orphan read mapped read in READ1! $prevname when $printcount lines are already printed, and $datalinecounter reads are red from the file.\n" ;
		$orphancountR1++ ;
	    }
	    
	    if ($waitingForR2==1) {
		print STDOUT "Discarding orphan read mapped read in READ2! $prevname when $printcount lines are already printed, and $datalinecounter reads are red from the file.\n" ;
		$orphancountR2++ ;
	    }
	    
	    $waitingForR2 = 0 ;
	    $waitingForR1 = 0 ;
	}
	
	#-----------------------------------------------------------------------------------------
	# Special cases, when we don't want to print but SKIP the read..
	#-----------------------------------------------------------------------------------------
	
	# If we already have both reads red..
	if ($flag & $READ1flag){
	    if ($pairprinted==1) { 
	    print STDOUT "Skipping multiple-mapped read in READ1! $name when $printcount lines are already printed, and $datalinecounter reads are red from the file.\n" ;
	    $duplicatecountR1++ ;
	    next LINE;
	    }
	}
	elsif ($flag & $READ2flag){
	    if ($pairprinted==1) { 
	    print STDOUT "Skipping multiple-mapped read in READ2! $name when $printcount lines are already printed, and $datalinecounter reads are red from the file.\n" ;
	    $duplicatecountR2++ ;
	    next LINE;
	    }
	}
	# If we have weird read not belonging to either R1 or R2
	else {
	    print STDOUT "Skipped a read not belonging in either of reads R1 or R2 ! $name when $printcount lines are already printed, and $datalinecounter reads are red from the file.\n" ;
	    next LINE;
	}
	
	#-----------------------------------------------------------------------------------------
	# Normal routine - if we will eventually print the read..
	#-----------------------------------------------------------------------------------------
	
	#Reverse complementing if needed
	
	if ($flag & $REV_COMPLflag){
	    #Reversing
	    my $seq = scalar reverse ("$thisline[9]") ;
	    my $score = scalar reverse ("$thisline[10]") ;
	    
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
	    
	    $seq =~ tr/ATCG/atcg/ ; #to be on the safe side, making this in 2 steps
	    $seq =~ tr/atcg/TAGC/ ;
	    
	    #Saving data
	    $data = "$seq\n+\n$score";    
	}
	else { # Saving data without modifications
	    my $seq = "$thisline[9]" ;
	    my $score = "$thisline[10]" ;
	    $data = "$seq\n+\n$score";    
	}
        
        if ($flag & $READ1flag){ # If the marker bit for READ1 is lit, then the data comes from READ1
	    
	    # We save the data to wait for printing
	    $dataR1inWaiting = "\@$name\n$data\n" ;
	    
	    # If R1 was what we were waiting for, we print.
	    if ($waitingForR1==1){
		&doPrinting ;
	    }
	    else{
		$waitingForR2=1;
	    }
	}
	
	if ($flag & $READ2flag){ # If the marker bit for READ2 is lit, then the data comes from READ2
	    
	    # We save the data to wait for printing
	    $dataR2inWaiting = "\@$name\n$data\n" ;
	    
	    # If R2 was what we were waiting for, we print.
	    if ($waitingForR2==1){
		&doPrinting ;
	    }
	    else{
		$waitingForR1=1;
	    }
	}
	
	
    }
}

# Update orphan counts if the last line was orphan..
	    if ($waitingForR1==1) {
		print STDOUT "Discarding orphan read mapped read in READ1! $prevname when $printcount lines are already printed, and $datalinecounter reads are red from the file.\n" ;
		$orphancountR1++ ;
	    }
	    
	    if ($waitingForR2==1) {
		print STDOUT "Discarding orphan read mapped read in READ2! $prevname when $printcount lines are already printed, and $datalinecounter reads are red from the file.\n" ;
		$orphancountR2++ ;
	    }




# In the end we wish to keep on eye what the code is doing ..
print STDOUT "TOTAL data entry lines read : $datalinecounter \n" ;
print STDOUT "Entries written to R1 and R2 output files (each got this amount of reads) : $printcount \n" ;

print STDOUT "Orphan reads found (and discarded before printing) : \n    Read1 had $orphancountR1 orphans, Read2 had $orphancountR2 orphans.\n" ;
print STDOUT "Duplicate-mapped reads found (and discarded before printing) : \n    Read1 had $duplicatecountR1 duplicates, Read2 had $duplicatecountR2 duplicates.\n" ;

#-----------------------------

# Subroutines

sub doPrinting{
    
 	    print PRINTER1 $dataR1inWaiting ;
   	    print PRINTER2 $dataR2inWaiting ;
	    
	    $waitingForR1 = 0;
	    $waitingForR2 = 0;
	    
    	    $pairprinted=1;
	    $prevname=$name;
	    $printcount++ ;
    
}