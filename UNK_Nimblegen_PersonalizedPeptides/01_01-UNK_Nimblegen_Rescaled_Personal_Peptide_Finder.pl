#!/usr/bin/perl

use strict;
use warnings;

####################################
# This script is to find personal peptides in UNK Nimblegen peptide array files.
####################################

my $line;
my @linearray;

use FileHandle;

# Process both rescaled and un-rescaled files.
my @files = ("Nimblegen_All_UNK_Log2_QN_lm_id_Rescaled.txt");
my @outfiles = ("Nimblegen_All_UNK_Log2_QN_lm_id_Rescaled_refpsn.txt", "Nimblegen_All_UNK_Log2_QN_lm_id_Rescaled_ref.txt", "Nimblegen_All_UNK_Log2_QN_lm_id_Rescaled_psn.txt");
open(IN, "$files[0]") || die;
open(OUT1, "> $outfiles[0]") || die;
open(OUT2, "> $outfiles[1]") || die;
open(OUT3, "> $outfiles[2]") || die;

my %peptides = (); # Hash containing all lines by protein names.
my %personalmark; # Hash to mark whether a peptide has personal counterpart.

# Print Header.
$line = <IN>;
print OUT1 $line;
print OUT2 $line;
print OUT3 $line;

# Now processing the file.
while($line = <IN>) {
	@linearray = split("\t", $line);
	my $key = $linearray[0];
	$personalmark{$key} = 0;
	my @linearray0 = split(";", $linearray[0]);
	if ($linearray0[0] =~ /:/) {
		my @linearray00 = split(":", $linearray0[0]);
		$key = $linearray00[0].";".$linearray0[1];
		$personalmark{$key} = 1;
	}
	push(@{$peptides{$key}}, $line);
}

# Now output the peptides with personal counterparts.
my @keys = sort {$a cmp $b} (keys %peptides);
foreach my $sortedkey (@keys) {
	next if ($personalmark{$sortedkey} != 1);
	if (scalar @{$peptides{$sortedkey}} == 2) {
		foreach my $finalrecord (@{$peptides{$sortedkey}}) {
			print OUT1 $finalrecord;
			@linearray = split("\t", $finalrecord);
			if ($linearray[0] !~ /:/) {
				print OUT2 $finalrecord;
			} else {
				print OUT3 $finalrecord;
			}
		}
	} else {
		my $markref = 0; # Marks whether this is the first element, which is ref peptide.
		my $firstref;
		foreach my $finalrecord (@{$peptides{$sortedkey}}) {
			@linearray = split("\t", $finalrecord);
			if ($markref == 0 && $linearray[0] !~ /:/) {$firstref = $finalrecord; $markref += 1;}
			if ($linearray[0] =~ /:/) {
				print OUT2 $firstref;
				print OUT3 $finalrecord;
				print OUT1 $firstref;
				print OUT1 $finalrecord;
			}
		}
	}
}

close IN;
close OUT1;
close OUT2;
close OUT3;

exit;
