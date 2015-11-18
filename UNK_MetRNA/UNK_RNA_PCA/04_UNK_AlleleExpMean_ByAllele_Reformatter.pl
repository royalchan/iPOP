#!/usr/bin/perl

use strict;
use warnings;

####################################
# This script is to reformat allele count mean by each allele in UNK_AlleleExpMean_NimblgenVar_MAF.txt.
# Ref information is based on the CADD v1.1 database.
# Created by Rui Chen, 11/03/2015.
####################################

my $line;
my @linearray;

# Read in Ref / Alt information from UNK variant file TEST-sny_vcf_RefAltInfo.txt.
open(VCF, "TEST-snp_vcf_RefAltInfo.txt") || die;
my %unkref;
my %unkalt;
my %hethom;
$line = <VCF>; # Get Rid of Header.
while ($line = <VCF>) {
	chomp $line;
	@linearray = split("\t", $line);
	my $key = $linearray[0]."_".$linearray[1];
	$unkref{$key} = $linearray[2];
	$unkalt{$key} = $linearray[3];
	if ($linearray[4] =~ /^0\/1/ || $linearray[4] =~ /^1\/0/) {
		$hethom{$key} = "Het";
	} elsif ($linearray[4] =~ /^1\/1/ || $linearray[4] =~ /^0\/0/) {
		$hethom{$key} = "Hom";
	} else {
		$hethom{$key} = "Unknown";
	}
}
close VCF;

# Now process the the real file.
open(IN, "UNK_AlleleExpMean_NimblgenVar_MAF.txt") || die;
open(OUT, "> UNK_AlleleExpMeanByAllele_NimblgenVar_MAF.txt") || die;

$line = <IN>; # Get Header.
print OUT "Chr\tPosition\tMaternal\/Paternal\tChrom_Pos_Ref_Alt\tESP_EUR_MAF\t1000G_EUR_MAF\tRef\tAllele_Count_Mean_Ref_hg19\tAllele_Count_Mean_Ref_mat\tAllele_Count_Mean_Ref_pat\tAllele_Count_Mean_Ref_mp\tAlt\tAllele_Count_Mean_Alt_hg19\tAllele_Count_Mean_Alt_mat\tAllele_Count_Mean_Alt_pat\tAllele_Count_Mean_Alt_mp\tHetHom\n";
while ($line = <IN>) {

	chomp $line;
	@linearray = split("\t", $line);
	my $key = $linearray[0]."_".$linearray[1]; # Chr_Coord.
	# Replace "." in MAF with 0 thus all are now numbers.
	if ($linearray[8] eq ".") {$linearray[8] = 0;}
	if ($linearray[9] eq ".") {$linearray[9] = 0;}
	
	# First determine Ref and Alt (Note that for chr12:108956433 ref/ alt disagrees between CAD and Mat/Pat).
	my @linearray7;
	my $ref;
	my $alt;
	if (defined $unkref{$key} || 0) {$ref = $unkref{$key};}
	if (defined $unkalt{$key} || 0) {$alt = $unkalt{$key};}
	if ((not defined $unkref{$key}) && ($linearray[7] !~ ".")) {
		@linearray7 = split("_", $linearray[7]);
		$ref = $linearray7[2];
		$alt = $linearray7[3];
	} elsif ((not defined $unkref{$key}) && ($linearray[7] =~ "."))  {
		$linearray[2] =~ s/Ref_//;
		$linearray[2] =~ s/Alt_//;
		@linearray7 =  split("\\/", $linearray[2]);
		$ref = $linearray7[0];
		$alt = $linearray7[1];
	}
	
	# Next perform reformatting.
	# Define an array with elements in the order consistent with header for output.
	my @output;
	$output[0] = $linearray[0]; # Chr.
	$output[1] = $linearray[1]; # Coord.
	$output[2] = $linearray[2]; # Maternal/Paternal.
	$output[3] = $linearray[7]; # Chrom_Pos_Ref_Alt.
	$output[4] = $linearray[8]; # ESP_EUR_MAF.
	$output[5] = $linearray[9]; # 1000G_EUR_MAF.
	# Het/Hom.
	if (defined $hethom{$key} || 0) {
		$output[16] = $hethom{$key};
	} else {
		$output[16] = "Unknown";
	}
	
	# Determine Ref / Alt.
	$linearray[2] =~ s/Ref_//;
	$linearray[2] =~ s/Alt_//;
	my @linearraygt = split("\\/", $linearray[2]);
	if ($linearraygt[0] eq $ref) {
		$output[6] = $linearraygt[0];
		$output[11] = $linearraygt[1];
	} else {
		$output[6] = $linearraygt[1];
		$output[11] = $linearraygt[0];
	}

	# hg19.
	my @linearraysub = split("\\:", $linearray[3]);
	my @linearrayct = split("\\/", $linearraysub[1]);
	if ($linearraygt[0] eq $ref) {
		$output[7] = $linearrayct[0];
		$output[12] = $linearrayct[1];
	} else {
		$output[7] = $linearrayct[1];
		$output[12] = $linearrayct[0];
	}

	# mat.
	@linearraysub = split("\\:", $linearray[4]);
	@linearrayct = split("\\/", $linearraysub[1]);
	if ($linearraygt[0] eq $ref) {
		$output[8] = $linearrayct[0];
		$output[13] = $linearrayct[1];
	} else {
		$output[8] = $linearrayct[1];
		$output[13] = $linearrayct[0];
	}

	# pat.
	@linearraysub = split("\\:", $linearray[5]);
	@linearrayct = split("\\/", $linearraysub[1]);
	if ($linearraygt[0] eq $ref) {
		$output[9] = $linearrayct[0];
		$output[14] = $linearrayct[1];
	} else {
		$output[9] = $linearrayct[1];
		$output[14] = $linearrayct[0];
	}

	# mp.
	@linearraysub = split("\\:", $linearray[6]);
	@linearrayct = split("\\/", $linearraysub[1]);
	if ($linearraygt[0] eq $ref) {
		$output[10] = $linearrayct[0];
		$output[15] = $linearrayct[1];
	} else {
		$output[10] = $linearrayct[1];
		$output[15] = $linearrayct[0];
	}

	# Print reformatted line.
	my $reline = join("\t", @output);
	print OUT "$reline\n";
}

close IN;
close OUT;

exit;
