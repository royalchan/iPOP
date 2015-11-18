#!/usr/bin/perl

use strict;
use warnings;

####################################
# This script is to find Minor Allele Frequencies for the personal variants in UNK Family in the CADD v1.1 file.
# Created by Rui Chen, 08/26/2015.
####################################

my $line;
my @linearray;

# Get information  from the CADD list I compiled.
open(CADD, "whole_genome_SNVs_inclAnno_AFEUR_slimmed_unq.tsv") || die;
# use Chr_Coord_Ref_Alt as key.
my %maf; # EUR MAF.
$line = <CADD>; # Get rid of Header.
$line = <CADD>; # Get rid of Header.
while ($line = <CADD>) {
	chomp $line;
	@linearray = split("\t", $line);
	my $key = $linearray[0]."_".$linearray[1];
	my $keylong = $linearray[0]."_".$linearray[1]."_".$linearray[2]."_".$linearray[3];
	my $espmaf = $linearray[4];
	if ($espmaf =~ /NA/) {$espmaf = 0;}
	my $tgmaf = $linearray[5];
	if ($tgmaf =~ /NA/) {$tgmaf = 0;}
	$maf{$key} = "$keylong\t$espmaf\t$tgmaf";
}
close CADD;

open(IN, "UNK_ALL_hg19matpatmp_norep_het_m_to_mnp_ratio_full.txt") || die;
open(OUT, "> UNK_ALL_hg19matpatmp_norep_het_m_to_mnp_ratio_full_MAF.txt") || die;
$line = <IN>; # Get Header.
chomp $line;
my $header = $line."\tChrom_Pos_Ref_Alt\tESP_EUR_MAF\t1000G_EUR_MAF\n";
print OUT $header;
while ($line = <IN>) {
	chomp $line;
	@linearray = split("\t", $line);
	my $key = $linearray[0]."_".$linearray[1];
	if (not defined $maf{$key} || 0) {$maf{$key} = ".\t.\t.";}
	print OUT "$line\t$maf{$key}\n";
}
close IN;
close OUT;

exit;
