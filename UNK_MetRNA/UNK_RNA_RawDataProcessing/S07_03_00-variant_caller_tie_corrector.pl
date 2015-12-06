#!/usr/bin/perl

# use strict;
# use warnings;

# I just found that my S07_02-variant_caller_from_mpileup.pl script mis-report the genotype of ties.
# This script is to mend the genotype and hom/het for the tie calls.
# I modified this script so it can be called by S07_03_00-LOADER_variant_caller_tie_corrector.pl and submitted
# 	paralleled jobs.
# Two command line arguments needed: $ARGV[0], input file name; $ARGV[1], output file name.
# Updated by Rui Chen, on September 13th, 2014.

use FileHandle;

my $dir = "/srv/gsfs0/projects/gbsc/iPoP/Users/danxie/unkrna/UALLMPILEUP/";

# Process file by file.
my $fh = FileHandle -> new;
my $oh = FileHandle -> new;

my $line;
my @linearray;

# Process the sorted files.
my $filename = $dir.$ARGV[0];
my $outname = $dir.$ARGV[1];
	
	$fh -> open("$filename");
	$oh -> open(">$outname");
	
	my $header = <$fh>;
	print $oh $header;
	
	while ($line = <$fh>) {
		
		chomp $line;
		@linearray =  split("\t", $line);
		
		my @counts = ($linearray[5], $linearray[6], $linearray[7], $linearray[8], $linearray[13], $linearray[14]);
		my @alleles = ("A", "T", "G", "C", $linearray[9], $linearray[10]);
		my @order = sort {$b <=> $a} @counts;
		
		# No need to change anything if line does not have ties.
		if ($order[0] != $order[1]) {
			print $oh "$line\n";
			next;
		}
		
		# Fix ties for the lines that have them.
		my @top2 = (); # Numbers for top 2 alleles.
		for (my $i = 0; $i < scalar (@counts); ++$i) {
			if ($counts[$i] == $order[0]) {
				push(@top2, $i);
			}
		}
		
		# Now re-determine the values for Het_Hom and Genotype.
		$linearray[2] = "het";
		$linearray[4] = "$alleles[$top2[0]]\/$alleles[$top2[1]]";
		
		# Finally let's rejoin the line and print it into the output file.
		my $reline = join("\t", @linearray);
		print $oh "$reline\n";
	
	}

	$fh -> close;
	$oh -> close;

exit;
