#!/usr/bin/perl

use strict;
use warnings;

####################################
# This script is to calculate allele count mean for each mapping from the UNK_AlleleExpression_NimblgenVar_MAF.txt file.
# Created by Rui Chen, 11/03/2015.
####################################

my $line;
my @linearray;

open(IN, "UNK_AlleleExpression_NimblgenVar_MAF.txt") || die;
open(OUT, "> UNK_AlleleExpMean_NimblgenVar_MAF.txt") || die;

$line = <IN>; # Get Header.
print OUT "Chr\tPosition\tMaternal\/Paternal\tAllele_Count_Mean_hg19\tAllele_Count_Mean_mat\tAllele_Count_Mean_pat\tAllele_Count_Mean_mp\tChrom_Pos_Ref_Alt\tESP_EUR_MAF\t1000G_EUR_MAF\n";
while ($line = <IN>) {

	@linearray = split("\t", $line);
	
	my $meanhg = ""; # hg19 mean in the format of M(Ref)/P(Alt):X/Y.
	my $meanmat = ""; # mat mean in the format of M(Ref)/P(Alt):X/Y.
	my $meanpat = ""; # pat mean in the format of M(Ref)/P(Alt):X/Y.
	my $meanmp = ""; # mp mean in the format of M(Ref)/P(Alt):X/Y.
	
	# Calculate hg19 mean.
	my $sum1 = 0; # Sum of first genotype (Mat or Ref).
	my $sum2 = 0; # Sum of second genotype (Pat or Alt).
	my $count = 0; # Sample number count.
	my $gt = ""; # Genotype.
	for (my $i = 67; $i <= 130; ++$i) {
		next if ($linearray[$i] =~ /NA/);
		++$count;
		my @linearraycell = split("\\:", $linearray[$i]);
		$gt = $linearraycell[0];
		my @linearraycount = split("\\/", $linearraycell[1]);
		$sum1 += $linearraycount[0];
		$sum2 += $linearraycount[1];
	}
	$sum1 /= $count;
	$sum2 /= $count;
	$meanhg = $gt.":".$sum1."/".$sum2;
	
	# Calculate mat mean.
	$sum1 = 0; # Sum of first genotype (Mat or Ref).
	$sum2 = 0; # Sum of second genotype (Pat or Alt).
	$count = 0; # Sample number count.
	$gt = ""; # Genotype.
	for (my $i = 198; $i <= 261; ++$i) {
		next if ($linearray[$i] =~ /NA/);
		++$count;
		my @linearraycell = split("\\:", $linearray[$i]);
		$gt = $linearraycell[0];
		my @linearraycount = split("\\/", $linearraycell[1]);
		$sum1 += $linearraycount[0];
		$sum2 += $linearraycount[1];
	}
	$sum1 /= $count;
	$sum2 /= $count;
	$meanmat = $gt.":".$sum1."/".$sum2;
	
	# Calculate pat mean.
	$sum1 = 0; # Sum of first genotype (Mat or Ref).
	$sum2 = 0; # Sum of second genotype (Pat or Alt).
	$count = 0; # Sample number count.
	$gt = ""; # Genotype.
	for (my $i = 329; $i <= 392; ++$i) {
		next if ($linearray[$i] =~ /NA/);
		++$count;
		my @linearraycell = split("\\:", $linearray[$i]);
		$gt = $linearraycell[0];
		my @linearraycount = split("\\/", $linearraycell[1]);
		$sum1 += $linearraycount[0];
		$sum2 += $linearraycount[1];
	}
	$sum1 /= $count;
	$sum2 /= $count;
	$meanpat = $gt.":".$sum1."/".$sum2;
	
	# Calculate mp mean.
	$sum1 = 0; # Sum of first genotype (Mat or Ref).
	$sum2 = 0; # Sum of second genotype (Pat or Alt).
	$count = 0; # Sample number count.
	$gt = ""; # Genotype.
	for (my $i = 460; $i <= 523; ++$i) {
		next if ($linearray[$i] =~ /NA/);
		++$count;
		my @linearraycell = split("\\:", $linearray[$i]);
		$gt = $linearraycell[0];
		my @linearraycount = split("\\/", $linearraycell[1]);
		$sum1 += $linearraycount[0];
		$sum2 += $linearraycount[1];
	}
	$sum1 /= $count;
	$sum2 /= $count;
	$meanmp = $gt.":".$sum1."/".$sum2;
	
	print OUT "$linearray[0]\t$linearray[1]\t$linearray[2]\t$meanhg\t$meanmat\t$meanpat\t$meanmp\t$linearray[524]\t$linearray[525]\t$linearray[526]";
}

close IN;
close OUT;

exit;
