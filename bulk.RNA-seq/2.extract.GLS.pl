#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <sid.dir> [sid.dir ...]\n\n";
	print STDERR "Extract the expression of GLS and KGA/GAC isoforms from RSEM results (supports human/mouse data).\n\n";
	exit 2;
}

my $minTPM = 1;	## only calculate ratios where GLS is higly expressed
my $fix = 0.1;

my $hs_KGA = 'NM_014905';
my $hs_GAC = 'NM_001256310';
my $mm_KGA = 'NM_001081081';
my $mm_GAC = 'NM_001113383';

my $hs_Q9UI32 = 'NM_013267';
my $mm_Q9UI32 = 'NM_001033264';

print "#Sid\tGLS\tKGA\tGAC\tLog2(KGA/GAC)\n";
## NOTE: KGA/GAC for aerobic glycolysis

foreach my $sid ( @ARGV ) {
	unless( -d $sid ) {
		print STDERR "ERROR: $sid is NOT a directory, SKIP!\n";
		next;
	}

	if( -s "$sid/rsem.genes.results" ) {
		open IN, "$sid/rsem.genes.results" or die( "$!" );
	} elsif( -s "$sid/$sid.rsem.genes.results" ) {
		open IN, "$sid/$sid.rsem.genes.results" or die( "$!" );
	} else {
		print STDERR "ERROR: $sid does not contain RSEM results!\n";
		next;
	}

	my ($GLS, $KGA, $GAC) = ( 0, 0, 0 );

	while( <IN> ) {
		chomp;
		my @l = split /\t/;	##gene_id	transcript_id(s)	length	effective_length	expected_count	TPM	FPKM
		if( $l[0] =~ /^GLS$/i ) {
			$GLS = $l[-2];	## use TPM
			last;
		}
	}
	close IN;

	if( -s "$sid/rsem.isoforms.results" ) {
		open IN, "$sid/rsem.isoforms.results" or die( "$!" );
	} else {
		open IN, "$sid/$sid.rsem.isoforms.results" or die( "$!" );
	}
	while( <IN> ) {
		chomp;
		my @l = split /\t/; ##transcript_id	gene_id	length	effective_length	expected_count	TPM	FPKM	IsoPct
		$l[0] =~ s/\..*//;
		if( $l[1] =~ /^GLS$/i ) {
			if( $l[0] eq $hs_KGA || $l[0] eq $mm_KGA ) {
				$KGA = $l[-3];
			} elsif( $l[0] eq $hs_GAC || $l[0] eq $mm_GAC ) {
				$GAC = $l[-3];
			}
		}
	}
	close IN;

	my $ratio = ($GLS > $minTPM) ? log(($KGA+$fix)/($GAC+$fix))/log(2) : 'NA';
	print join("\t", $sid, $GLS, $KGA, $GAC, $ratio), "\n";
}

