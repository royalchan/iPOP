#!/usr/bin/perl

use strict;
use warnings;

# This script is to generate assumed cytoscape correlation network input
#	from the UNK_NRSOC_Spearman_RhoSquared.txt file.

open(IN, "UNK_RNASVA_Spearman_RhoSquared.txt") || die "Cannot open the input file!!\n\n";
open(OUT, "> UNK_RNASVA_Spearman_RhoSquared_CytoscapeInput.txt");

my $line;
my @linearray;

$line = <IN>; # Get Time Point Name Column.
chomp $line;
my @colnames = split("\t", $line);
shift(@colnames);
print OUT "UNK_TP_1\tUNK_TP_2\tSpearman_Rho_Squared\n";

my $i = 0;
while ($line = <IN>) {
	chomp $line;
	@linearray = split("\t", $line);
	my $rowname = shift(@linearray);
	for (my $j = $i + 1; $j < scalar @linearray; ++$j) {
		if ($linearray[$j] >= 0.16) {
			print OUT "$rowname\t$colnames[$j]\t$linearray[$j]\n";
		}
	}
	++$i;
}

close IN;
close OUT;

open(IN, "UNK_RNASVA_Spearman_Rho.txt") || die "Cannot open the input file!!\n\n";
open(OUT, "> UNK_RNASVA_Spearman_Rho_CytoscapeInput.txt");

$line = <IN>; # Get Time Point Name Column.
chomp $line;
@colnames = split("\t", $line);
shift(@colnames);
print OUT "UNK_TP_1\tUNK_TP_2\tSpearman_Rho\n";

$i = 0;
while ($line = <IN>) {
	chomp $line;
	@linearray = split("\t", $line);
	my $rowname = shift(@linearray);
	for (my $j = $i + 1; $j < scalar @linearray; ++$j) {
		if ($linearray[$j] >= 0.4 || $linearray[$j] <= -0.4) {
			print OUT "$rowname\t$colnames[$j]\t$linearray[$j]\n";
		}
	}
	++$i;
}

close IN;
close OUT;

exit;
