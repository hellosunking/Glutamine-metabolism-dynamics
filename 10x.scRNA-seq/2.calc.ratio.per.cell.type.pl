#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 1 ) {
	print STDERR "\nUsage: $0 <GLS1.stat> <cell.anno> [target.column=1] [cell.prefix=null]\n";
	print STDERR "\nNote: columns start from 0; if you have more than 1 column, e.g., 1:2, 1:2:3\n\n";
	exit 2;
}

my $target= $ARGV[2] || 1;
my @column = split /:/, $target;

my $cellprefix = $ARGV[3] || '';

my $fix = 1;	## fix in log2-normalization
my $log2 = log(2);

my (%GAC, %KGA);
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;	#Barcode GAC KGA
	$l[0] =~ s/-1$//;
	$l[0] = "$cellprefix$l[0]";	## add prefix, useful after pooling
	$GAC{$l[0]} = $l[1];
	$KGA{$l[0]} = $l[2];
}
close IN;

my (%all, %valid, %G, %K);

open IN, "$ARGV[1]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/; #Barcode Anno1 Anno2 ...
	next if $#l < $column[-1];

	$l[0] =~ s/-1$//;
	my @annotation = map { $l[$_] } @column;
	my $anno = join(":", @annotation);
	$all{$anno} ++;
	if( exists $GAC{$l[0]} ) {
		$valid{$anno} ++;
		$G{$anno} += $GAC{$l[0]};
		$K{$anno} += $KGA{$l[0]};
	}
}
close IN;

print "#Type\tAll.ID\tValid.ID\tGAC\tKGA\tLog2(KGA/GAC)\n";
foreach my $c ( sort keys %all ) {
	my $GAC = $G{$c} || 0;
	my $KGA = $K{$c} || 0;
	print join("\t", $c, $all{$c}, $valid{$c}||0, $GAC, $KGA, log(($KGA+$fix)/($GAC+$fix))/$log2), "\n";
}

