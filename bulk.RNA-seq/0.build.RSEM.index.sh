#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 3 ]
then
	echo "Usage: $0 <genes.gtf> <genome.fa> <index.id>"
	exit 2
fi > /dev/stderr

RefSeq=$1
fasta=$2
genome=$3

thread=`cat /proc/cpuinfo | grep processor | wc -l`
rsem-prepare-reference -p $thread --bowtie2 --gtf $RefSeq $fasta $genome

