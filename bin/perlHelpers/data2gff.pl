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
#use warnings ;

require Time::localtime ;
my $timestamp ;

use File::Basename;
use Pod::Usage;
use Getopt::Long;

my $counter = 0 ;
my $datalinecounter = 0 ;
my $printcount = 0 ;

my $noheading = 0 ;
my $input = "";
my $inputName = "";
my $source = "";
my $type = "";
my $output = "";
my $test = "";

my $columnString = "";
my $colnameString = "";

my @columns ;
my @colnames ;

my $sourceType = "";
my $col_6to8 = ".\t.\t.";
my @facetNames = "";

my @thisline ;

my $headingRed = 0 ;

{
#Creating unique names for output files
my @time_array = localtime(time) ;
$time_array[5] = $time_array[5]-100 ;
$time_array[4] = $time_array[4]+1 ;
$timestamp = "$time_array[3]$time_array[4]$time_array[5]_$time_array[2]$time_array[1]$time_array[0]"
}

#Automated help message when the program was called funny
my $message = <<END;
\n
Reads in any tab-delimited file (3 first columns Chr Start Stop), with or without heading. Transforms it to .gff file. Comment lines in input file should begin with #.
\n
perl data2gff.pl -data <InputFile.fileExtension>
    [-test <only print out names and column numbers of selected data> ]
    [-source <experiment name etc> to be printed to .gff f.ex \"MyCellLine\"]
    [-type <subtype for the data> to be printed to .gff f.ex \"MacsData\"]
    [--col <columns to be printed to FACETS in OutPut.gff file> f.ex: -col 1,2 prints 4th and 5th column of input file - first 3 columns are Chr Start Stop which are always printed. Default : ALL DATA] 
    [--colnames <names for columns to be printed to FACETS in OutPut.gff file> f.ex -colnames p_value,pileup  Default : use names provided in input file] 
    [-noheading <No column header line in .txt file> Default : FALSE = use names provided in input file>]
    [-output <output.gff> Default : InputFileName_TimeStamp.gff]

Does not check for heading line formatting - expects tab-delimited fields, one per column, no special characters. If non-standard heading, prints "" instead of the name of the column - then the column names can be set with -colnames option.
    
END

pod2usage("$0: $message" ) if (($#ARGV<0) && (-t STDIN));

#Reading ARGV in with a built-in method in Getopt package
&GetOptions ("data=s"    => \$input,
	     "source:s"	 => \$source,
	     "type:s"	 => \$type,
             "test"    => \$test,
             "col:s"     => \$columnString,
             "colnames:s"=> \$colnameString,
             "noheading" => \$noheading,
	     "output:s"  => \$output
	  );

# Printing version name
print "";
print "Perl version :";
print "$^V";
print "";

#Setting input-name
$inputName = basename($input) ;

#Setting output name, if it was not given in arguments
if ($output){
    $output = ">".$output ;
}
else {
    my $name = $inputName;    
    $name =~ s/\./_/;
    $output = ">$name"."_"."$timestamp".".gff";
}

#Writing down which colums are to be printed out
if ($columnString){
    @columns = split(/,/, $columnString) ;
#    print "\$columnString = $columnString\n\@columns = @columns\n"
}

#Writing down the names of the columns to be printed
if ($colnameString){
    @colnames = split(/,/,$colnameString);
    
#    print "\$colnameString = $colnameString\n\@colnames = @colnames\n"
}


#Formatting columns 2 and 3 (data source, data type) for printing (they will be the same for each data line)
#Also removing non-aplhanumeric characters from the names - they may mix up the loading of the .gff file to SQL database

if ($source){
    $source =~ s/!\w/_/g ;
}
if ($type){
    $type =~ s/!\w/_/g ;
}


if($source && $type){
    $sourceType = "$source\t$type" ;
}
elsif ($source){
    $sourceType = "$source\t." ;
}
elsif ($type){
    $sourceType = ".\t$type" ;
}
else {
    $sourceType = ".\t." ;
}



open(DATA, $input) || die $!;

#Opening output - unless we are only testing (via -test option), in which case no output is written to a file
unless ($test){
open(PRINTER, $output) || die $!; 
print STDOUT "Beginning to read data..\n" ;
}

LINE : while (<DATA>) {
    
    @thisline = "";
    
    chomp;
    
    # Heading line is red different - if the column names are to be red from the heading line, it is allowed to be either space or tab delimited
#    if (!$noheading && !$headingRed ) {
#	@thisline=split/\s/ ;
 #   }
    # Data lines are tab-delimited
  #  else {
	@thisline=split/\t/ ;
   # }
    
    $counter++ ;
    if ( $counter % 1000 == 0){
    print STDOUT "One thousand lines processed\n";    
    }
    

    {
    #If we have data line, not a comment line (or empty line or just whitespace)
    unless (($thisline[0] =~ /^#/) || ($thisline[0] eq "") || ($thisline[0] =~ /^\s$/)) {
        
        #Heading is captured here
        if ($datalinecounter==0){
	    
	  
	    
	    
	    #If we print ONLY CERTAIN COLUMNS
	    if ($columnString){
		#If user has given names for these columns in command line
		if ($colnameString){
		    foreach my $Col (0..$#columns){
		    ($facetNames[$Col] = $colnames[$Col]) or die "ARGUMENT ERROR : There should be equally many FACETs mentioned in -col and -colname !" ;
		    }
		}
		#If no column names available neither from command line nor from input file
		elsif ($noheading){
		    foreach my $Col(0..$#columns){
		    $facetNames[$Col] = "data_". ($Col+1) ;
		    }
		}
		#Reading column names from input file (default)
		else {
		  
		    foreach my $Col (0..$#columns){
		    $facetNames[$Col] = $thisline[$columns[$Col]+2] or die "ERROR : There is not a proper header line in the input file !";
		    

		    
		    }
		}
		
	    }
	    
	    #If we print ALL COLUMNS (default)
	    else {
		#If user has given names for these columns in command line
		if ($colnameString){
		    foreach my $Col (0..$#colnames){
		    ($facetNames[$Col] = $colnames[$Col])or die "ARGUMENT ERROR : There should be equally many FACETs mentioned in -colname as there are data fields (in addition to \"chr\" \"start\" and \"stop\") in InputFile.xls!" ;
		    }    
		}
		#If no column names available neither from command line nor from input file
		elsif ($noheading){
		    foreach my $Col(3..$#thisline){
		    $facetNames[$Col-3] = "data_". ($Col-2) ;
		    }
		}
		#Reading column names from input file (default)
		else {
		    foreach my $Col (3..$#thisline){
		    $facetNames[$Col-3] = $thisline[$Col] or die "ERROR : There is not a proper header line in the input file !"; 
		    }    
		}
		
	    }
	    
	    #Removing non-aplhanumeric characters from FACET names - they may mix up the loading of the .gff file to SQL database
	    
	    foreach my $Col (0..$#facetNames) {
	    $facetNames[$Col] =~ s/\W/_/g ; # replaces everything non-alphanumeric with _
	    }
	    
	    $headingRed = 1 ;
	    
	    
	    
	    
	    
	    #Printing out what is to be outputted in FACETs :
	    
	    if(@columns){
		print STDOUT "FACETs to be outputted [InputColNo, FacetName]:\n" ;
		foreach my $colNo (0..$#columns){
		    my $printColNo = $columns[$colNo] ;
		    print STDOUT "[ $printColNo , $facetNames[$colNo] ]\n" ;
		}
	    }
	    else {
		print STDOUT "FACETs to be outputted [FacetNo, FacetName]:\n" ;
		foreach my $facetNo (0..$#facetNames) {
		    my $printColNo = $facetNo + 1 ;
		    print STDOUT "[ $printColNo , $facetNames[$facetNo] ]\n";
		}
	    }
	    
	    #If the outputting of printed-to-be FACETs was all we wanted - i.e. if "-test" option is on :
	    last LINE if ($test);

	    #If we are reading the column amounts in a DATA LINE, we will also save this data.
	    #Otherwise jump to the beginning of the loop, and start reading real data
	    unless ($noheading){
	    $datalinecounter++;
	    next LINE ;
	    }
        }
	
        
	
        # Printing neatly - preparing data !
	#
	# macs2.xls format : chr, start, end, length, abs_summit, pileup, -LOG10(pvalue), fold_enrichment, -LOG10(qvalue), name
	#
	# GFF3 FORMAT :	chromosome, source, type, start, stop, xx, +/-strand, xx, FACETS
	# Our format :	chromosome, $source_type, start, stop, $col_6to8, FACETS

	
	my $col_1to8 ;
	{
	#my $chr   = @thisline[0];
	#my $start = @thisline[1];
	#my $stop  = @thisline[2];
	
	my $chr   = shift @thisline;
	my $start = shift @thisline;
	my $stop  = shift @thisline;
	
	$col_1to8 = "$chr\t$sourceType\t$start\t$stop\t$col_6to8" ;
	}
	
	my $facets = "";
	
	#If only certain columns are to be printed
	if (@columns){
	    foreach my $facetNo(0..$#columns){
		my $name  = "$facetNames[$facetNo]" ;
		my $value = "$thisline[$columns[$facetNo]-1]" ;
		$facets = "$facets"."$name"."="."$value"."; "
	    }
	}	    
	
	#If all columns are to be printed
	else{
	    foreach my $facetNo(0..$#facetNames){
		my $name  = "$facetNames[$facetNo]";
		my $value = "$thisline[$facetNo]" ;
		$facets = "$facets"."$name"."="."$value"."; "
	    }
	}
	
	#Printing the data
        
	select PRINTER ;
	print PRINTER "$col_1to8\t$facets\n";
	
	
	$datalinecounter++;
	next LINE ;
	
	}
    }
}

close DATA ;

#Printing out the summary of the running of the script. Unless it was a test (and thus no output created) - then nothing is printed.
unless ($test){
close PRINTER ;
my $dataLinesMinusHeader ;
if ($noheading){
    $dataLinesMinusHeader = $datalinecounter ;
}
else {
   $dataLinesMinusHeader = $datalinecounter -1 ; 
}
print STDOUT "Total amount lines red : $counter \nTotal data lines (including heading line) red : $datalinecounter and printed (gff does not have header line): $dataLinesMinusHeader \n" ;
}
