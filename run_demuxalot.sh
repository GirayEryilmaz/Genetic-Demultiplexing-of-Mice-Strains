#!/usr/bin/env bash

#SBATCH --job-name=demuxalot
#SBATCH --output=demuxalot_%j.log
#SBATCH --partition=compute
#SBATCH --qos=batch
#SBATCH --time=26:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH	--cpus-per-task=70
#SBATCH --mem=200GB

SAMPLE='MYSAMPLE'

OUTDIR="$SAMPLE/demuxalot"

BARCODES="$SAMPLE/cellranger/filtered_feature_bc_matrix/barcodes.tsv.gz"

BAM="$SAMPLE/cellranger/possorted_genome_bam.bam"
INDS="strains.txt"

VCF="B6_NZO_CAST.fully.reordered.vcf.gz"


python run_Demuxalot.py \
        -b $BARCODES \
        -a $BAM \
        -n $INDS \
        -v $VCF \
        -o $OUTDIR \
        -p 70 \
        -r True