#!/usr/bin/perl

# use strict;
# use warnings;

# This script attempts to call diploid Variant/Indel-containing genotype from the .mpileup file generated with SAMTools.
# This script is to be called by the S07_02-LOADER_variant_caller_from_mpileup.pl file.
# Two variables are passed as $ARGV: $ARGV[0] -- Input file name; $ARGV[1] -- Output file name.
# Created by Rui Chen, on July 28th, 2014.
# Updated by Rui Chen, on August 1st, 2014.
# 	Updates include: A) Added exclusion criteria for output:
#						A1) Only count reads with mapping quality >= 15 (ord(QUALITY_ASCII) -33 >= 15); 
#						A2)	Homozygous: Allele read count >=20; Heterozygous: Minor allele read count >= 10;
# 					 B) Changed the order of output columns, and added genotype and mean quality score to output.
#					 C) Added an Error output file "UALL.mpileup.count.err".

use FileHandle;
use List::Util qw(sum);

my $dir = "/srv/gsfs0/projects/gbsc/iPoP/Users/danxie/unkrna/UALLMPILEUP/";

my $filename = $dir.$ARGV[0];
my $outname = $dir.$ARGV[1];
my $errname = $dir."UALL.mpileup.count.err";

# Process file by file.
my $fh = FileHandle -> new;
my $oh = FileHandle -> new;

my $line;
my @linearray;

$fh -> open("$filename");
$oh -> open(">$outname");
open (ERR, ">>$errname");
	
# my $header = "Chr\tPosition\tHet_Hom\tRef\tGenotype\tA_Count\tT_Count\tG_Count\tC_Count\tIndel_1\tIndel_2\tIndel_1_Length\tIndel_2_Length\tIndel_1_Count\tIndel_2_Count\tA_AverageBaseQuality\tT_AverageBaseQuality\tG_AverageBaseQuality\tC_AverageBaseQuality\tIndel_1_AverageBaseQuality\tIndel_2_AverageBaseQuality\n";
my $header = "Chr\tPosition\tHet_Hom\tRef\tGenotype\tA_Count\tT_Count\tG_Count\tC_Count\tIndel_1\tIndel_2\tIndel_1_Length\tIndel_2_Length\tIndel_1_Count\tIndel_2_Count\tA_AverageBaseQuality\tT_AverageBaseQuality\tG_AverageBaseQuality\tC_AverageBaseQuality\n";
print $oh $header;
	
while ($line = <$fh>) {
	chomp $line;
	@linearray =  split("\t", $line);
		
#	next if ($linearray[4] !~ /[\+\-]\d+/);

#	Write lines with scalar(@linearray) != 6 directly to ERR.
	if ($linearray[3] != 0 && scalar(@linearray) != 6) {
		print ERR "For file $ARGV[0], column number is not 6 at line:\n$line\n\n";
		next;
	}

	# Quality Score Penalty.
	my @qualpenalty = ();
		
	my %indels; # Hash of indels in the active row, with the indel genotype as key, and count number as value.
	my %indelslength;
	# The following variable is for the actual Indel documentation: [\+\-]\d+.{\d+}
	my $indeltoadd;
	my $linearray4 = $linearray[4]; # Base Info.
	my $linearray5 = $linearray[5]; # Qual Info.
	
	# Count Indels first.
	# Assume each position could have at most 2 different indels from maternal and paternal alleles, 
	#	therefore they will be replaced by "X" (First Identified) and "Y" (Second Identified).
	my @indelgenotype1; # Documentation of the identified indel genotypes in the order by which they are identified.
	my $k = 0; # Counter of @indelgenotype1.
	my $linearray4uc = uc $linearray4;
	while ($linearray4uc =~ /[\+\-]\d+/) {
		$indeltoadd = $linearray4uc;
		my $stringlength = $indeltoadd;
		my $signedstringlength = $indeltoadd;
		$stringlength =~ s/^.*[\+\-](\d+)[ATGCatgcNn]+.*$/$1/;
		$signedstringlength =~ s/^.*([\+\-]\d+)[ATGCatgcNn]+.*$/$1/;
		$indeltoadd =~ s/^.*([\+\-]\d+\D{$stringlength}).*$/$1/;
		$indelgenotype1[$k] = $indeltoadd;
		my $replacement;
		if ($indeltoadd eq $indelgenotype1[0]) {
			$replacement = "X";
		} else {$replacement = "Y";}
		$linearray4uc =~ s/\Q$indeltoadd\E/\Q$replacement\E/;
# This step was cancelled as it could not distinguish opposite indels such as +2CC and -2CC, though rare.
# 		# Further restrict $indeltoadd to indel bases.
# 		$indeltoadd =~ s/^.*[\+\-]\d+(\D{$stringlength}).*$/$1/;
# 		$indeltoadd = uc $indeltoadd;
		$indels{$indeltoadd} += 1;
		$indelslength{$indeltoadd} = $signedstringlength;
		$indeltoadd = "";
		++$k;
	}
	my $firstindelgenytype = $indelgenotype1[0]; # Genotype of first indel without length and indel information.
	$firstindelgenytype = uc $firstindelgenytype;
	my @indelcount = sort {$b <=> $a} (values %indels); # Sorted counts for each indel in this row.
	my $keymax = ""; # key for the max value for %indels.
	my $keysecond = ""; # key for the second largest value for %indels.
	my $keymaxshort = ""; # short form of key for the max value for %indels.
	my $keysecondshort = ""; # short form of key for the second largest value for %indels.
	my $indelcount2 = 0; # Read count for Potential 2nd Indel, this variable is to deal especially with the case when there is no $indelcount[1].
	if (scalar @indelcount == 1) {
		foreach (keys %indels) {
				$keymax = $_;
				$keymaxshort = $keymax;
				my $stringlength = $keymaxshort;
				$stringlength =~ s/^[\+\-](\d+)[ATGCatgcNn]+$/$1/;
				$keymaxshort =~ s/^[\+\-]\d+(\D{$stringlength})$/$1/;
				$keymaxshort = uc $keymaxshort;
				$indelcount2 = 0;
				$keysecond = $keysecondshort = "--";
				$indelslength{$keysecond} = 0;
		}
	} elsif (scalar @indelcount > 1) {
		foreach (keys %indels) {
			if ($indels{$_} eq $indelcount[0]) {
				$keymax = $_;
				$keymaxshort = $keymax;
				my $stringlength = $keymaxshort;
				$stringlength =~ s/^[\+\-](\d+)[ATGCatgcNn]+$/$1/;
				$keymaxshort =~ s/^[\+\-]\d+(\D{$stringlength})$/$1/;
				$keymaxshort = uc $keymaxshort;
			}
			if ($indels{$_} eq $indelcount[1]) {
				$keysecond = $_;
				$keysecondshort = $keysecond;
				my $stringlength = $keysecondshort;
				$stringlength =~ s/^[\+\-](\d+)[ATGCatgcNn]+$/$1/;
				$keysecondshort =~ s/^[\+\-]\d+(\D{$stringlength})$/$1/;
				$keysecondshort = uc $keysecondshort;
			}
		$indelcount2 = $indelcount[1];
		}
	} else {
		$keymax = $keymaxshort = "--";
		$keysecond = $keysecondshort = "--";
		$indelslength{$keymax} = 0;
		$indelslength{$keysecond} = 0;
		$indelcount[0] = 0;
		$indelcount2 = 0;
	}
			
	# Count of ATGC.
	my $linearray4foradj = $linearray4uc; # For base quality adjustment.
	my $refcount = ($linearray4foradj =~ tr/[\,\.]//);
	my $acount = ($linearray4foradj =~ tr/[Aa]//);
	my $tcount = ($linearray4foradj =~ tr/[Tt]//);
	my $gcount = ($linearray4foradj =~ tr/[Gg]//);
	my $ccount = ($linearray4foradj =~ tr/[Cc]//);
	if ($linearray[2] =~ /A/) {
		$acount = $refcount;
	} elsif ($linearray[2] =~ /T/) {
		$tcount = $refcount;
	} elsif ($linearray[2] =~ /G/) {
		$gcount = $refcount;
	} elsif ($linearray[2] =~ /C/) {
		$ccount = $refcount;
	}
	
	# Count Adjustment by base quality.
	# NOTE: INDELS DO NOT HAVE BASE QUALITY!!!
	# Further cleaning up of $linearray4foradj.
	$linearray4foradj =~ s/\^.//g;
	$linearray4foradj =~ s/\$//g;
	$linearray4foradj =~ s/\./\Q$linearray[2]\E/g;
	$linearray4foradj =~ s/\,/\Q$linearray[2]\E/g;	
	my @bases = split("", $linearray4foradj);
	# BIG DISCOVERY!! INDELS DO NOT HAVE CORRESPONDING BASE QUALITY SCORES!!!
	# Correct $linearray5 with "!" (zero quality score) by positions of indels ("X" or "Y").
	my @positions; # Positions of indels ("X" or "Y").
	for (my $i = 0; $i < scalar @bases; ++$i) {
		if ($bases[$i] eq "X" || $bases[$i] eq "Y") {
			push(@positions, $i);
		}
	}
	my @readqualori = split("", $linearray5);
	my @readqual; # Base quality array post correction.
	if (scalar @positions == 0) {
		@readqual = @readqualori;
	} else {
		my $j = 0; # Independent Counter of Real (base) Positions.
		my $firstposition = shift(@positions); # Position shifted from @positions.
		for (my $i = 0; $i < scalar @readqualori; ++$i) {
			while ($j eq $firstposition) {
				push(@readqual, "!");
				++$j;
				$firstposition = shift(@positions);
			}
			push(@readqual, $readqualori[$i]);
			++$j;
		}
		# Deal with the last bases (in case they are a string of indels).
		# Fix the situation when X and/or Y are at the end of @bases.
		while ($firstposition) {
			if ($j eq $firstposition) {
					push(@readqual, "!");
					++$j;
			} elsif ($j < $firstposition) {
					do {
						push(@readqual, "!");
						++$j;
					} until ($j eq $firstposition);
			}
			$firstposition = shift(@positions);
		}
	}

	# Deduct any count with base quality < 15.
	my @phred = basequaltonumber(@readqual);
	if ((scalar @phred) != (scalar @bases)) {
		my $qn = scalar @phred;
		my $bn = scalar @bases;
		print ERR "For file $ARGV[0], Base number ($bn) and Quality number ($qn) DO NOT Match at line:\n$line\n\n";
	}
# 	my @indel1phred;
# 	my $indel1failedcount = 0;
# 	my @indel2phred;
# 	my $indel2failedcount = 0;
	my @aphred;
	my $afailedcount = 0;
	my @tphred;
	my $tfailedcount = 0;
	my @gphred;
	my $gfailedcount = 0;
	my @cphred;
	my $cfailedcount = 0;
	for (my $i = 0; $i < scalar(@phred); ++$i) {
		if ($phred[$i] < 15) {
# 			if ($bases[$i] eq "X") {
# 				++$indel1failedcount;
# 			} elsif ($bases[$i] eq "Y") {
# 				++$indel2failedcount;
# 			} elsif ($bases[$i] eq "A") {
			if ($bases[$i] eq "A") {
				++$afailedcount;
			} elsif ($bases[$i] eq "T") {
				++$tfailedcount;
			} elsif ($bases[$i] eq "G") {
				++$gfailedcount;
			} elsif ($bases[$i] eq "C") {
				++$cfailedcount;
			}
		} else {
# 			if ($bases[$i] eq "X") {
# 				$indel1phred[$i] = $phred[$i];
# 			} elsif ($bases[$i] eq "Y") {
# 				$indel2phred[$i] = $phred[$i];
# 			} elsif ($bases[$i] eq "A") {
			if ($bases[$i] eq "A") {
				push(@aphred, $phred[$i]);
			} elsif ($bases[$i] eq "T") {
				push(@tphred, $phred[$i]);
			} elsif ($bases[$i] eq "G") {
				push(@gphred, $phred[$i]);
			} elsif ($bases[$i] eq "C") {
				push(@cphred, $phred[$i]);
			}		
		}
	}
	# Correct base counts by subtracting failed counts.
# 	if (scalar @indelcount >= 1) {
# 		if ($keymax eq $firstindelgenytype) {
# 			$indelcount[0] = $indelcount[0] - $indel1failedcount;
# 			$indelcount2 = $indelcount2 - $indel2failedcount;
# 		} elsif ($keysecond eq $firstindelgenytype) {
# 			$indelcount[0] = $indelcount[0] - $indel2failedcount;
# 			$indelcount2 = $indelcount2 - $indel1failedcount;
# 		}
# 	}
	$acount = $acount - $afailedcount;
	$tcount = $tcount - $tfailedcount;
	$gcount = $gcount - $gfailedcount;
	$ccount = $ccount - $cfailedcount;
# 	if ($indelcount[0] < 0) {$indelcount[0] = 0;}
# 	if ($indelcount2 < 0) {$indelcount2 = 0;}
	if ($acount < 0) {$acount = 0;}
	if ($tcount < 0) {$tcount = 0;}
	if ($gcount < 0) {$gcount = 0;}
	if ($ccount < 0) {$ccount = 0;}
	# Calculate Average Base Quality Score.
# 	my $aveindelcount1maq = 0;
# 	my $aveindelcount2maq = 0;
# 	if (scalar @indelcount == 1) {
# 		if (scalar @indel1phred > 0) {
# 			$aveindelcount1maq = sum(@indel1phred) / @indel1phred;
# 		}
# 		$aveindelcount2maq = 0;
# 	}
# 	if (scalar @indelcount >= 2) {
# 		if ($keymax eq $firstindelgenytype) {
# 			if (scalar @indel1phred > 0) {
# 				$aveindelcount1maq = sum(@indel1phred) / @indel1phred;
# 			}
# 			if (scalar @indel2phred > 0) {
# 				$aveindelcount2maq = sum(@indel2phred) / @indel2phred;
# 			}
# 		} elsif ($keysecond eq $firstindelgenytype) {
# 			if (scalar @indel2phred > 0) {
# 				$aveindelcount1maq = sum(@indel2phred) / @indel2phred;
# 			}
# 			if (scalar @indel1phred > 0) {
# 				$aveindelcount2maq = sum(@indel1phred) / @indel1phred;
# 			}
# 		}
# 	}
	my $aveamaq = 0;
	my $avetmaq = 0;
	my $avegmaq = 0;
	my $avecmaq = 0;	
	if (scalar @aphred > 0) {
		$aveamaq = sum(@aphred) / @aphred;
	}
	if (scalar @tphred > 0) {
		$avetmaq = sum(@tphred) / @tphred;
	}
	if (scalar @gphred > 0) {
		$avegmaq = sum(@gphred) / @gphred;
	}
	if (scalar @cphred > 0) {
		$avecmaq = sum(@cphred) / @cphred;
	}
	
	# Print to output.
	# Do not output if allele has a read depth < 20 for homozygous position, or minor allele has a read depth < 10 for heterozygous position.
	my @genotypeorderedbyreadcountsort = sort {$b <=> $a} (($indelcount[0], $indelcount2, $acount, $tcount, $gcount, $ccount)); # Order genotype by read depth from high to low.
	my %match = (
		$indelcount[0] => $keymax,
		$indelcount2 => $keysecond,
		$acount => "A",
		$tcount => "T",
		$gcount => "G",
		$ccount => "C",
	);
	next if ($genotypeorderedbyreadcountsort[0] < 20 && $genotypeorderedbyreadcountsort[1] < 10);
	my $hethom; # Value for het or hom calls.
	my $genotype; # Called Genotype of the position.
	if ($genotypeorderedbyreadcountsort[1] == 0) {
		$hethom = "hom";
		$genotype = "$match{$genotypeorderedbyreadcountsort[0]}/$match{$genotypeorderedbyreadcountsort[0]}";
	} elsif ($genotypeorderedbyreadcountsort[1] > 0 && $genotypeorderedbyreadcountsort[0] / $genotypeorderedbyreadcountsort[1] > 2) {
		$hethom = "hom";
		$genotype = "$match{$genotypeorderedbyreadcountsort[0]}/$match{$genotypeorderedbyreadcountsort[0]}";
	} else {
		$hethom = "het";
		$genotype = "$match{$genotypeorderedbyreadcountsort[0]}/$match{$genotypeorderedbyreadcountsort[1]}";
	} 
	
#	print $oh "$linearray[0]\t$linearray[1]\t$hethom\t$linearray[2]\t$genotype\t$acount\t$tcount\t$gcount\t$ccount\t$keymaxshort\t$keysecondshort\t$indelslength{$keymax}\t$indelslength{$keysecond}\t$indelcount[0]\t$indelcount2\t$aveamaq\t$avetmaq\t$avegmaq\t$avecmaq\t$aveindelcount1maq\t$aveindelcount2maq\n";
	print $oh "$linearray[0]\t$linearray[1]\t$hethom\t$linearray[2]\t$genotype\t$acount\t$tcount\t$gcount\t$ccount\t$keymaxshort\t$keysecondshort\t$indelslength{$keymax}\t$indelslength{$keysecond}\t$indelcount[0]\t$indelcount2\t$aveamaq\t$avetmaq\t$avegmaq\t$avecmaq\n";
}

$fh -> close;
$oh -> close;
close ERR;

exit;

############
# SUBROUTINE
############
# Convert ASCII Base Quality Scores to Numeric Phred Scores.
sub basequaltonumber {

	my @basequal = @_; # Input Base Quality Score.
	my @numbers; # Output Phred Score in Numbers.
	foreach (@basequal) {
		my $num = ord($_) - 33;
		push(@numbers, $num);
	}
	return @numbers;

}