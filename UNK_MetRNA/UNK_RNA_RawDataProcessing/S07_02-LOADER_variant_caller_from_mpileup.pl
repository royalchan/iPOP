#!/usr/bin/perl

# use strict;
# use warnings;

# This script attempts to load the S07_02-variant_caller_from_mpileup.pl file in submitted jobs,
#	and call diploid Indel-containing genotype from the .mpileup file generated with SAMTools.
# Created by Rui Chen, on July 28th, 2014.

my @filename; # Name of the input file (mapped with hg19).
my @filenamem; # Name of the input file (maternal).
my @filenamep; # Name of the input file (paternal).
my @filenamerd; # Name of the input file (mapped with hg19, duplicate reads removed).
my @filenamerdm; # Name of the input file (maternal, duplicate reads removed).
my @filenamerdp; # Name of the input file (paternal, duplicate reads removed).
my @outname; # Name of the output file (mapped with hg19).
my @outnamem; # Name of the output file (maternal).
my @outnamep; # Name of the output file (paternal).
my @outnamerd; # Name of the output file (mapped with hg19, duplicate reads removed).
my @outnamerdm; # Name of the output file (maternal, duplicate reads removed).
my @outnamerdp; # Name of the output file (paternal, duplicate reads removed).

# Read in all the input files and determine the output file names.
# Read in file list.
open (FILELIST, "S07-variantcallerfilelist.txt") or die "Cannot open the file list!!!\n\n";
my $i = 0;
while (my $line = <FILELIST>) {
	chomp $line;
	my $noreprmdup; # Determine whether to use norep or rmdup in the file names.
	if ($line =~ /U42RNA/ || $line =~ /U46RNA/ || $line =~ /U47RNA/ || $line =~ /U48RNA/ || $line =~ /U49RNA/ || $line =~ /U50RNA/ || $line =~ /U52RNA/ || $line =~ /U56RNA/) {
		$noreprmdup = "rmdup";
	} else {$noreprmdup = "norep";}
	$filename[$i] = $line.".sorted.bam.mpileup";
	$filenamem[$i] = $line."_m.sorted.bam.mpileup";
	$filenamep[$i] = $line."_p.sorted.bam.mpileup";
	$filenamerd[$i] = $line.".".$noreprmdup.".sorted.bam.mpileup";
	$filenamerdm[$i] = $line."_m.".$noreprmdup.".sorted.bam.mpileup";
	$filenamerdp[$i] = $line."_p.".$noreprmdup.".sorted.bam.mpileup";
	$outname[$i] = $line.".sorted.bam.mpileup.count";
	$outnamem[$i] = $line."_m.sorted.bam.mpileup.count";
	$outnamep[$i] = $line."_p.sorted.bam.mpileup.count";
	$outnamerd[$i] = $line.".".$noreprmdup.".sorted.bam.mpileup.count";
	$outnamerdm[$i] = $line."_m.".$noreprmdup.".sorted.bam.mpileup.count";
	$outnamerdp[$i] = $line."_p.".$noreprmdup.".sorted.bam.mpileup.count";
	++$i;
}

my $line;
my @linearray;

my $dir = "/srv/gsfs0/projects/gbsc/iPoP/Users/danxie/unkrna/UALLMPILEUP/";
my $dir2 = "/srv/gsfs0/projects/gbsc/iPoP/Users/danxie/unkrna/UALLMPILEUP/qeo/";

for (my $i = 0; $i < (scalar @filename); ++$i) {
	
	open (F1, "> $dir2$filename[$i]readcount.q");
	open (F2, ">> $dir2$filename[$i]readcount.e");
	open (F3, ">> $dir2$filename[$i]readcount.o");
	print F1 "#!/bin/tcsh\n\ndate\n";
	my $time = localtime;
	print F2 "\n############\n$time\n############\n";
	print F3 "\n############\n$time\n############\n";

	print F1 "perl $dir"."S07_02-variant_caller_from_mpileup.pl $filename[$i] $outname[$i]\n";

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
	my $commandII = "$dir2$filename[$i]readcount.o -e ";
	my $commandIII = "$dir2$filename[$i]readcount.e ";
	my $commandIV = "$dir2$filename[$i]readcount.q";
	my $command = $commandI.$commandII.$commandIII.$commandIV;
	print $command."\n";
	system("$command");

}

for (my $i = 0; $i < (scalar @filenamem); ++$i) {
	
	open (F1, "> $dir2$filenamem[$i]readcount.q");
	open (F2, ">> $dir2$filenamem[$i]readcount.e");
	open (F3, ">> $dir2$filenamem[$i]readcount.o");
	print F1 "#!/bin/tcsh\n\ndate\n";
	my $time = localtime;
	print F2 "\n############\n$time\n############\n";
	print F3 "\n############\n$time\n############\n";

	print F1 "perl $dir"."S07_02-variant_caller_from_mpileup.pl $filenamem[$i] $outnamem[$i]\n";

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
                   " -q extended -R y -N RC_RCm$i -o ";
	my $commandII = "$dir2$filenamem[$i]readcount.o -e ";
	my $commandIII = "$dir2$filenamem[$i]readcount.e ";
	my $commandIV = "$dir2$filenamem[$i]readcount.q";
	my $command = $commandI.$commandII.$commandIII.$commandIV;
	print $command."\n";
	system("$command");
	
}

for (my $i = 0; $i < (scalar @filenamep); ++$i) {
	
	open (F1, "> $dir2$filenamep[$i]readcount.q");
	open (F2, ">> $dir2$filenamep[$i]readcount.e");
	open (F3, ">> $dir2$filenamep[$i]readcount.o");
	print F1 "#!/bin/tcsh\n\ndate\n";
	my $time = localtime;
	print F2 "\n############\n$time\n############\n";
	print F3 "\n############\n$time\n############\n";

	print F1 "perl $dir"."S07_02-variant_caller_from_mpileup.pl $filenamep[$i] $outnamep[$i]\n";

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
                   " -q extended -R y -N RC_RCp$i -o ";
	my $commandII = "$dir2$filenamep[$i]readcount.o -e ";
	my $commandIII = "$dir2$filenamep[$i]readcount.e ";
	my $commandIV = "$dir2$filenamep[$i]readcount.q";
	my $command = $commandI.$commandII.$commandIII.$commandIV;
	print $command."\n";
	system("$command");
	
}

for (my $i = 0; $i < (scalar @filenamerd); ++$i) {
	
	open (F1, "> $dir2$filenamerd[$i]readcount.q");
	open (F2, ">> $dir2$filenamerd[$i]readcount.e");
	open (F3, ">> $dir2$filenamerd[$i]readcount.o");
	print F1 "#!/bin/tcsh\n\ndate\n";
	my $time = localtime;
	print F2 "\n############\n$time\n############\n";
	print F3 "\n############\n$time\n############\n";

	print F1 "perl $dir"."S07_02-variant_caller_from_mpileup.pl $filenamerd[$i] $outnamerd[$i]\n";

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
                   " -q extended -R y -N RC_RCrd$i -o ";
	my $commandII = "$dir2$filenamerd[$i]readcount.o -e ";
	my $commandIII = "$dir2$filenamerd[$i]readcount.e ";
	my $commandIV = "$dir2$filenamerd[$i]readcount.q";
	my $command = $commandI.$commandII.$commandIII.$commandIV;
	print $command."\n";
	system("$command");
	
}

for (my $i = 0; $i < (scalar @filenamerdm); ++$i) {
	
	open (F1, "> $dir2$filenamerdm[$i]readcount.q");
	open (F2, ">> $dir2$filenamerdm[$i]readcount.e");
	open (F3, ">> $dir2$filenamerdm[$i]readcount.o");
	print F1 "#!/bin/tcsh\n\ndate\n";
	my $time = localtime;
	print F2 "\n############\n$time\n############\n";
	print F3 "\n############\n$time\n############\n";

	print F1 "perl $dir"."S07_02-variant_caller_from_mpileup.pl $filenamerdm[$i] $outnamerdm[$i]\n";

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
                   " -q extended -R y -N RC_RCrdm$i -o ";
	my $commandII = "$dir2$filenamerdm[$i]readcount.o -e ";
	my $commandIII = "$dir2$filenamerdm[$i]readcount.e ";
	my $commandIV = "$dir2$filenamerdm[$i]readcount.q";
	my $command = $commandI.$commandII.$commandIII.$commandIV;
	print $command."\n";
	system("$command");
	
}

for (my $i = 0; $i < (scalar @filenamerdp); ++$i) {
	
	open (F1, "> $dir2$filenamerdp[$i]readcount.q");
	open (F2, ">> $dir2$filenamerdp[$i]readcount.e");
	open (F3, ">> $dir2$filenamerdp[$i]readcount.o");
	print F1 "#!/bin/tcsh\n\ndate\n";
	my $time = localtime;
	print F2 "\n############\n$time\n############\n";
	print F3 "\n############\n$time\n############\n";

	print F1 "perl $dir"."S07_02-variant_caller_from_mpileup.pl $filenamerdp[$i] $outnamerdp[$i]\n";

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
                   " -q extended -R y -N RC_RCrdp$i -o ";
	my $commandII = "$dir2$filenamerdp[$i]readcount.o -e ";
	my $commandIII = "$dir2$filenamerdp[$i]readcount.e ";
	my $commandIV = "$dir2$filenamerdp[$i]readcount.q";
	my $command = $commandI.$commandII.$commandIII.$commandIV;
	print $command."\n";
	system("$command");
	
}

exit;
