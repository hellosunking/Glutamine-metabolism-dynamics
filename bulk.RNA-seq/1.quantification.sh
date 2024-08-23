#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 1 ]
then
	echo "
Usage: $0 [options] -i <RSEM.index.prefix> -1 <R1.fq> [-2 <R2.fq>] -o <output.sample.id>

This pipeline is designed for handling bulk RNA-seq or SMART-seq data.

Options:
  -t thread  Set running threads. Default: all threads (>=4 required)
  -k kit     Set kit for trimming adaptors. Default: illumina. Supports bgi, illumina
  -L info    Specify the strandness of the library. Supports 'unstranded' (default), 'stranded' or 'forward', and 'reverse'.

Source codes of Ktrim and Krmdup could be found at:
https://github.com/hellosunking/Ktrim/
https://github.com/hellosunking/Microcket/
"
	exit 2
fi > /dev/stderr

currSHELL=`readlink -f $0`
PRG=`dirname $currSHELL`

## default parameters
read1=""
read2=""
sid=""

RSEMindex=""
thread=0
seqKit="illumina"
strandess="unstranded"
rsem_strand=""

## read command line parameters
while getopts ":i:o:1:2:t:k:L" OPTION
do
	case $OPTION in
		i)RSEMindex="$OPTARG"
			;;
		o)sid="$OPTARG"
			;;
		1)read1="$OPTARG"
			;;
		2)read2="$OPTARG"
			;;
		t)thread="$OPTARG"
			;;
		k)seqKit="$OPTARG"
			;;
		L)strandess="$OPTARG"
			;;
		?)echo -e "\n\n***** ERROR: unsupported option detected. *****\n"
			;;
	esac
done

## check parameters
if [ -z "$read1" ] || [ -z "$sid" ] || [ -z "$RSEMindex" ]
then
	echo "Error: compulsory parameters are missing!"
	exit 101
fi

## check other parameters
if [ $thread == 0 ]     ## use all threads
then
	thread=`cat /proc/cpuinfo | grep processor | wc -l`
elif [ $thread -lt 4 ]
then
	echo "ERROR: at least 4 threads are needed."
	exit 11
fi

if [ $strandess == "unstranded" ]
then
	rsem_strand="-forward-prob 0.5"
elif [ $strandess == "stranded" ] || [ $strandess == "forward" ]
then
	rsem_strand="-forward-prob 1"
elif [ $strandess == "reverse" ] || [ $strandess == "reversed" ]
then
	rsem_strand="-forward-prob 0"
else
	echo "ERROR: unsupported strandness! Must be unstranded/stranded/reverse!"
	exit 12
fi

## start analysis
echo "Preprocessing using Kombo ..."
if [ -z "$read2" ]	## single-end
then
	$PRG/ktrim -k $seqKit -U $read1 -o $sid -t $thread -c | $PRG/krmdup.se -i /dev/stdin -o $sid.rmdup

	echo "Running RSEM ..."
	mkdir -p $sid
	cd $sid
	rsem-calculate-expression $rsem_strand -p $thread --bowtie2 --no-bam-output \
		../$sid.rmdup.read1.fq $RSEMindex $sid.rsem
	cd ../
	rm -f $sid.rmdup.read1.fq
else	## paired-end
	$PRG/ktrim -k $seqKit -t $thread -o $sid -1 $read1 -2 $read2 -c | $PRG/krmdup.pe -i /dev/stdin -o $sid.rmdup

	echo "Running RSEM ..."
	mkdir -p $sid
	cd $sid
	rsem-calculate-expression $rsem_strand -p $thread --bowtie2 --no-bam-output \
		--paired-end ../$sid.rmdup.read1.fq ../$sid.rmdup.read2.fq $RSEMindex $sid.rsem
	cd ../
	rm -f $sid.rmdup.read1.fq $sid.rmdup.read2.fq
fi

mv $sid.rmdup.log $sid.trim.log $sid
echo "Done. The outputs are stored in '$sid' directory."

