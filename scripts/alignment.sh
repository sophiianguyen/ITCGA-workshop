#!/bin/bash
#SBATCH --job-name=Alignment # you can give your job a name
#SBATCH --nodes 1 # the number of processors or tasks
#SBATCH --cpus-per-task=4
#SBATCH --account=itcga # our account
#SBATCH --time=10:00:00 # the maximum time for the job
#SBATCH --mem=20gb # the amount of RAM
#SBATCH --partition=itcga # the specific server in chimera we are using
#SBATCH --reservation=2025_JUNE_ITCGA_WORKSHOP #This reservation lets us cut in line to use itcga cores
#SBATCH --error=%x-%A_%a.err   # a filename to save error messages into
#SBATCH --output=%x-%A_%a.out  # a filename to save any printed output into

#input_dir=$1 # takes this from the command line, first item after the script
#output_dir=$2 # takes this from the command line, second item



module load gcc-10.2.0-gcc-9.3.0-f3oaqv7
module load python-3.8.12-gcc-10.2.0-oe4tgov
module load hisat2-2.1.0-gcc-9.3.0-u7zbyow
module load samtools-1.10-gcc-9.3.0-flukja5
module load subread-2.0.2-gcc-10.2.0

for i in ~/project/data/trim_fastq/KAS_fastq/*.lite.1.fastq
do

name=$(basename ${i} .lite.1.fastq)

hisat2 -x ~/1_project/data/genome/hg38/hg38 \
-U ~/project/data/trim_fastq/KAS_fastq/${name}.lite.1.fastq \
-S ~/project/results/sam/${name}.sam \
-p 4

samtools view -S -b ~/project/results/sam/${name}.sam > ~/project/results/bam/${name}.bam
samtools sort -o ~/project/results/bam/${name}_sorted.bam ~/project/results/bam/${name}.bam
samtools flagstat ~/project/results/bam/${name}.bam
samtools index ~/project/results/bam/${name}_sorted.bam


done
featureCounts -a ~/1_project/data/genome/Homo_sapiens.GRCh38.111.gtf \
-o ~/project/results/counts/counts.txt \
-T 24 \
~/project/results/bam/*_sorted.bam
