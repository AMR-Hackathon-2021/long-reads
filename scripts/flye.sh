#!/bin/bash

# Get options
OPT_SIZE=10
OPT_THREADS=8
OPT_OUTDIR=assemblies
while getopts ":o:s:t:d:" opt; do
    case $opt in
        o)
            OPT_OUTDIR=$OPTARG
            ;;
        s)
            OPT_SIZE=$OPTARG
            ;;
        t)
            OPT_THREADS=$OPTARG
            ;;
  
        d)
            debug=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done
INPUT_FILE=$1
BASENAME=$(basename $INPUT_FILE | cut -f1 d_ | cut -f1 -d.)

if [ ! -e "$INPUT_FILE"]; then
    echo "Input file $INPUT_FILE does not exist"
    exit 1
fi

set -euxo pipefail

mkdir -p $OPT_OUTDIR
flye --genome-size ${OPT_SIZE}M --nano-raw "$INPUT_FILE" --threads $OPT_THREADS -o $OPT_OUTDIR/$BASENAME
minimap2 -t $OPT_THREADS -ax map-ont $OPT_OUTDIR/$BASENAME/assembly.fasta "$INPUT_FILE" | samtools view -bS | samtools sort -@ 4  -o $OPT_OUTDIR/$BASENAME/aln.bam -
covtobed $OPT_OUTDIR/$BASENAME/aln.bam > $OPT_OUTDIR/$BASENAME/cov.bed
abricate $OPT_OUTDIR/$BASENAME/assembly.fasta > OPT_OUTDIR/$BASENAME/abricate.tab
