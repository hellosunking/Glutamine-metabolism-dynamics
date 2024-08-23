#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <possorted_genome_bam.bam> [species=hs|mm] [cell.prefix=NULL]\n";
	print STDERR "\nExtract the GLS isoform usage from cell-ranger output BAM.\n\n";
	exit 2;
}

my $minQ = 30;

my ($loci, $ENST_GAC, $ENST_KGA);
my $genome = $ARGV[1] || 'hs';
my $prefix = $ARGV[2] || "";

## GLS1 gene loci and isoform definition
if( $genome eq 'hs' ) {	## ENSG00000115419
	$loci = 'chr2:190,880,797-190,965,552';
	$ENST_GAC = 'ENST00000338435';
	$ENST_KGA = 'ENST00000320717';
} elsif( $genome eq 'mm' ) {	## ENSMUSG00000026103
	$loci = 'chr1:52,162,887-52,191,895';
	$ENST_GAC = 'ENSMUST00000114510';
	$ENST_KGA = 'ENSMUST00000114513';
} else {
	print STDERR "ERROR: unknown species!\n";
	exit 1;
}

my %GAC;	##the shorter one (tumor preferred)
my %KGA;	##the longer one

open IN, "samtools view -@ 2 $ARGV[0] $loci |" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;
	next if ($l[4] < $minQ) || ($l[1] & 256) || ($l[1] & 512) || ($l[1] & 1024);	## low quality, or non-primary alignment result, QC failed, etc

	next unless /CR:Z:(\S+)/;
	my $barcode=$1;

	next unless /UR:Z:(\S+)/;
	my $umi=$1;

	if( /TX:Z:$ENST_GAC/ ) {
		next if /TX:Z:$ENST_KGA/;	## conflict

		$GAC{$barcode}->{"$umi:$l[3]"} = 1;
#		print STDERR "GAC: $barcode $umi:$l[3]\n";
	} elsif( /TX:Z:$ENST_KGA/ ) {
		$KGA{$barcode}->{"$umi:$l[3]"} = 1;
#		print STDERR "KGA: $barcode $umi:$l[3]\n";
	}
}
close IN;

print "#Barcode\tGAC\tKGA\n";
foreach my $barcode ( keys %GAC ) {
	my $c1 = keys %{$GAC{$barcode}};
	if( exists $KGA{$barcode} ) {
		my $c2 = keys %{$KGA{$barcode}};
		print "$prefix$barcode\t$c1\t$c2\n";
		delete $KGA{$barcode};
	} else {
		print "$prefix$barcode\t$c1\t0\n";
	}
}

foreach my $barcode ( keys %KGA ) {
	my $c2 = keys %{$KGA{$barcode}};
	print "$prefix$barcode\t0\t$c2\n";
}

