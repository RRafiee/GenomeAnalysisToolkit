#!/bin/bash -eu
#$ -cwd -V
#$ -pe smp 1
#$ -l h_rt=24:00:00
#$ -l h_vmem=6G
#$ -R y
#$ -q all.q,bigmem.q


# 2017
# Converts aligned BAM to gziped FASTQ using BAM file passed in at command-line.
# Uses SAMtools to randomly order reads in aligned BAM file, to avoid mapping
# bias as discussed here:
# http://gatkforums.broadinstitute.org/discussion/2908/howto-revert-a-bam-file-to-fastq-format

set -o pipefail
hostname
date

#source ../GATKsettings.sh

## Base dir - should auto set to where this script resides
BASE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

TMPDIR=$BASE_DIR
SAMTOOLS="/opt/software/bsu/bin/samtools"


BAM=$1
B_NAME=$(basename $BAM)

echo "** Variables **"
echo " - BASE_DIR = $BASE_DIR"
echo " - PWD = $PWD"
echo " - BAM = $BAM"
echo " - B_NAME = $B_NAME"

#read -p "before copy process"
#echo "Copying input $BAM to $TMPDIR/"
#/usr/bin/time --verbose cp -v $BAM $TMPDIR

echo "Running SAMtools on $TMPDIR/$B_NAME saving output as a gziped FASTQ file"
/usr/bin/time --verbose $SAMTOOLS bamshuf -uOn 128 $TMPDIR/$B_NAME $TMPDIR/tmp | $SAMTOOLS bam2fq - | gzip > $TMPDIR/$B_NAME.fastq.gz

echo "Copying $TMPDIR/*.fastq.gz to $PWD"
/usr/bin/time --verbose cp -v $TMPDIR/*.fastq.gz $PWD

echo "Deleting $TMPDIR/*.bam"
rm $TMPDIR/*.bam

echo "Deleting $TMPDIR/*.gz"
rm $TMPDIR/*.gz

date
echo "END"
