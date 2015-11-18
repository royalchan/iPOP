#!/usr/bin/perl

use strict;
use warnings;

####################################
# This script is to slim the UNK_ALL_hg19matpatmp_norep_het_m_to_mnp_ratio_full_MAF.txt file down to only the variants on Nimblegen Arrays.
# Created by Rui Chen, 11/03/2015.
####################################

my $line;
my @linearray;

# Get variant information from the UNK_Family_Genotype_MAF_ProbeID.txt file.
open(NIM, "UNK_Family_Genotype_MAF_ProbeID.txt") || die;
# use Chr_Coord as key.
my %coord; # EUR MAF.
$line = <NIM>; # Get rid of Header.
while ($line = <NIM>) {
	@linearray = split("\t", $line);
	my $key = $linearray[2]."_".$linearray[3];
	$key =~ s/^chr//;
	$coord{$key} = $key;
}
close NIM;

open(IN, "UNK_ALL_hg19matpatmp_norep_het_m_to_mnp_ratio_full_MAF.txt") || die;
open(OUT, "> UNK_AlleleExpression_NimblgenVar_MAF.txt") || die;
$line = <IN>; # Get Header.
print OUT $line;
while ($line = <IN>) {
	@linearray = split("\t", $line);
	my $key = $linearray[0]."_".$linearray[1];
	next if (not defined $coord{$key} || 0);
	print OUT $line;
}
close IN;
close OUT;

exit;
