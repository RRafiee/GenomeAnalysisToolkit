# Genome Analysis Toolkit Key Points

===============================================================================
# What is paired-end read in DNA sequencing?

Illumina paired-end sequencing is based on the idea that you have initial DNA fragments (longer than your actual read length) and you sequence both its ends. On the Illumina chip, both ends of each sequence are amplified prior to actual sequencing using bridging. This approach results in two reads per fragment, with the first read in forward orientation and the second read in reverse-complement orientation. Depending on the inital fragment size and read length, these fragment can either overlap or not. In Single-end reads, the sequence fragment are sequenced from one direction only.

FASTQ is a text file format (human readable) that provides 4 lines of data per sequence:
Sequence identifier,
The sequence,
Comments,
Quality scores.

FASTQ format is commonly used to store sequencing reads, in particular from Illumina and Ion Torrent platforms.
Paired-end reads may be stored either in one FASTQ file (alternating) or in two different FASTQ files. Paired-end reads may have sequence identifiers ended by "/1" and "/2" respectively. 
You can check the orientation of a FASTQ file by 'cat' command line in linux.
$ cat filename | more  (use CTRL+z to stop) 

===============================================================================
# How to revert a BAM file back to FastQ?

Revert a BAM file back to FastQ. This comes in handy when you receive data that has been processed but not according to GATK Best Practices, and you want to reset and reprocess it properly.

Prerequisites
Installed HTSlib
Steps
Shuffle the reads in the bam file
Revert the BAM file to FastQ format
Compress the FastQ file

## 1. Shuffle the reads in the bam file

Shuffle the reads in the bam file so they are not in a biased order before alignment by running the following HTSlib command:
htscmd bamshuf -uOn 128 aln_reads.bam tmp > shuffled_reads.bam 
This creates a new BAM file containing the original reads, which still retain their mapping information, but now they are no longer sorted. The aligner uses blocks of paired reads to estimate the insert size. If you don’t shuffle your original bam, the blocks of insert size will not be randomly distributed across the genome, rather they will all come from the same region, biasing the insert size calculation. This is a very important step which is unfortunately often overlooked.

## 2. Revert the BAM file to FastQ

Revert the BAM file to FastQ format by running the following HTSlib command:
htscmd bam2fq -a shuffled_reads.bam > interleaved_reads.fq 
This creates an interleaved FastQ file called interleaved_reads.fq containing the now-unmapped paired reads.
Interleaved simply means that for each pair of reads in your paired-end data set, both the forward and the reverse reads are in the same file, as opposed to having them in separate files.


## 3. Compress the FastQ file

Compress the FastQ file to reduce its size using the gzip utility:

gzip interleaved_reads.fq

This creates a gzipped FastQ file called interleaved_reads.fq.gz. This file is ready to be used as input for the Best Practices workflow. BWA handles gzipped fastq files natively, so you don’t need to unzip the file to use it later on.



