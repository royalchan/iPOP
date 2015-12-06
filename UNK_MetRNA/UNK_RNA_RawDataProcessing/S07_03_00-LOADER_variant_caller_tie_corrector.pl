#!/usr/bin/perl

# use strict;
# use warnings;

# This script is to load the S07_03_00-variant_caller_tie_corrector.pl file in submitted jobs,
#	and correct genotypes in lines with ties from the .mpileup file generated with SAMTools.
# Created by Rui Chen, on September 13th, 2014.

use FileHandle;

my $dir = "/srv/gsfs0/projects/gbsc/iPoP/Users/danxie/unkrna/UALLMPILEUP/";
my $dir2 = "/srv/gsfs0/projects/gbsc/iPoP/Users/danxie/unkrna/UALLMPILEUP/qeo/";

my $line;
my @linearray;

open(FILELIST, "S07-variantcallerfilelist.txt");
my @filename;
my @filenamend;
my @outname;
my @outnamend;
my $i = 0;
while ($line = <FILELIST>) {
	chomp $line;
	my $noreprmdup;
	if ($line =~ /U42RNA/ || $line =~ /U46RNA/ || $line =~ /U47RNA/ || $line =~ /U48RNA/ || $line =~ /U49RNA/ || $line =~ /U50RNA/ || $line =~ /U52RNA/ || $line =~ /U56RNA/) {
		$noreprmdup = "rmdup";
	} else {$noreprmdup = "norep";}
	$filename[$i] = $line.".sorted.bam.mpileup.count";
	$filenamend[$i] = $line.".".$noreprmdup.".sorted.bam.mpileup.count";
	$outname[$i] = $line.".sorted.bam.mpileup.count.tiecor";
	$outnamend[$i] = $line.".".$noreprmdup.".sorted.bam.mpileup.count.tiecor";
	++$i;
	$filename[$i] = $line."_m.sorted.bam.mpileup.count";
	$filenamend[$i] = $line."_m.".$noreprmdup.".sorted.bam.mpileup.count";
	$outname[$i] = $line."_m.sorted.bam.mpileup.count.tiecor";
	$outnamend[$i] = $line."_m.".$noreprmdup.".sorted.bam.mpileup.count.tiecor";
	++$i;
	$filename[$i] = $line."_p.sorted.bam.mpileup.count";
	$filenamend[$i] = $line."_p.".$noreprmdup.".sorted.bam.mpileup.count";
	$outname[$i] = $line."_p.sorted.bam.mpileup.count.tiecor";
	$outnamend[$i] = $line."_p.".$noreprmdup.".sorted.bam.mpileup.count.tiecor";
	++$i;
}
close FILELIST;
	
# Process the sorted files.
for (my $i = 0; $i < scalar (@filename); ++$i) {
		
	open (F1, "> $dir2$filename[$i]tiecor.q");
	open (F2, ">> $dir2$filename[$i]tiecor.e");
	open (F3, ">> $dir2$filename[$i]tiecor.o");
	print F1 "#!/bin/tcsh\n\ndate\n";
	my $time = localtime;
	print F2 "\n############\n$time\n############\n";
	print F3 "\n############\n$time\n############\n";

	print F1 "perl $dir"."S07_03_00-variant_caller_tie_corrector.pl $filename[$i] $outname[$i]\n";

	print F1 "date\necho fin.\n\n";

	close(F1);
	close(F2);
	close(F3);

	my $commandI = "qsub ".
	               "-m ea ".
                   "-l h_stack=10M -l h_vmem=3.5G ".
                   "-pe shm 1 -w e ".
           #       "-l h_vmem=13G ".
           #       "-pe shm 1 -w e ".
                   " -q extended -R y -N RC_RC$i -o ";
	my $commandII = "$dir2$filename[$i]tiecor.o -e ";
	my $commandIII = "$dir2$filename[$i]tiecor.e ";
	my $commandIV = "$dir2$filename[$i]tiecor.q";
	my $command = $commandI.$commandII.$commandIII.$commandIV;
	print $command."\n";
	system("$command");

}

# Process the norep sorted files.
for (my $i = 0; $i < scalar (@filenamend); ++$i) {
		
	open (F1, "> $dir2$filenamend[$i]tiecor.q");
	open (F2, ">> $dir2$filenamend[$i]tiecor.e");
	open (F3, ">> $dir2$filenamend[$i]tiecor.o");
	print F1 "#!/bin/tcsh\n\ndate\n";
	my $time = localtime;
	print F2 "\n############\n$time\n############\n";
	print F3 "\n############\n$time\n############\n";

	print F1 "perl $dir"."S07_03_00-variant_caller_tie_corrector.pl $filenamend[$i] $outnamend[$i]\n";

	print F1 "date\necho fin.\n\n";

	close(F1);
	close(F2);
	close(F3);

	my $commandI = "qsub ".
	               "-m ea ".
                   "-l h_stack=10M -l h_vmem=3.5G ".
                   "-pe shm 1 -w e ".
           #       "-l h_vmem=13G ".
           #       "-pe shm 1 -w e ".
                   " -q extended -R y -N RC_RC$i -o ";
	my $commandII = "$dir2$filenamend[$i]tiecor.o -e ";
	my $commandIII = "$dir2$filenamend[$i]tiecor.e ";
	my $commandIV = "$dir2$filenamend[$i]tiecor.q";
	my $command = $commandI.$commandII.$commandIII.$commandIV;
	print $command."\n";
	system("$command");

}

exit;
