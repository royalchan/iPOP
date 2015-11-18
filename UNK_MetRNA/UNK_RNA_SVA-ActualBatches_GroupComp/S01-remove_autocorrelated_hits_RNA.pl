#!/usr/bin/perl

use strict;
use warnings;

####################################
# Remove autocorrelated hits from ShapiroWilk filtered analytes.
####################################

my $line;
my @linearray;

# Index autocorrelated hits.
open(ACF, "UNK_RNASeq_hg19_QNSVA_LSPMOD_ALL_Original_ACF_Bootstrap_Hits.txt") || die;
my %acf;
$line = <ACF>; # Get rid of header.
while ($line = <ACF>) {
	@linearray = split("\t", $line);
	$acf{$linearray[0]} = 1;
}
close ACF;

# Stream original list and output non-autocorrelated hits.
open(IN, "UNK_RNASeq_NearRealBatch_ComBat_Rescaled.txt") || die;
open(OUT, "> UNK_RNASeq_NearRealBatch_ComBat_Rescaled_no_ac.txt") || die;
$line = <IN>;
print OUT $line; # Output header.
while ($line = <IN>) {
	@linearray = split("\t", $line);
	if (not defined $acf{$linearray[0]} || 0) {
		print OUT $line;
	}
}
close IN;
close OUT;

exit;
