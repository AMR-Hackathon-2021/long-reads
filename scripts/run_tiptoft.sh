#!/bin/bash
set -euo pipefail

# specify the data directory
DATADIR=/data/long-reads/data/MinION

# create output directory if not exists
OUTDIR=/data/long-reads/data/tiptoft

if [ ! -d "${OUTDIR}" ]; then
 mkdir -p $OUTDIR
fi

# run the job
for file in $(ls $DATADIR/*.fastq.gz); do
    sample_id=`basename ${file} | cut -f1 -d.`
    echo -e "predicting plasmid in $sample_id"
    tiptoft "$file" --output_file $OUTDIR/$sample_id.tiptoft.out
done
