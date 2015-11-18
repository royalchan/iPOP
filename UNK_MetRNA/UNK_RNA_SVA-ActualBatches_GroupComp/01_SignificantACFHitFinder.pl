#!/usr/bin/perl

use strict;
use warnings;

use FileHandle;

my $line;
my @linearray;

my $cutoff = 0.1819333;

opendir(CURR, ".") || die;
my @files = grep(/\.txt$/, readdir(CURR));
closedir CURR;
my @outfiles;
foreach (@files) {
	my $tmp = $_;
	$tmp =~ s/\.txt//;
	$tmp .= "_Hits.txt";
	push(@outfiles, $tmp);
}

my $fh = FileHandle -> new;
my $oh = FileHandle -> new;

for (my $i = 0; $i < scalar @files; ++$i) {

	$fh -> open("< $files[$i]");
	$oh -> open("> $outfiles[$i]");

	# Output header.
	$line = <$fh>;
	print $oh $line;
	
	while ($line = <$fh>) {
		@linearray = split("\t", $line);
		if ($linearray[1] > $cutoff) {print $oh $line;}
	}
	
	$fh -> close;
	$oh -> close;
	
}

exit;
