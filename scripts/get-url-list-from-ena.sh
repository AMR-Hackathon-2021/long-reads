#!/bin/bash

#https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB47011&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true&limit=0

URL_PREFIX="https://www.ebi.ac.uk/ena/portal/api/filereport?accession="
URL_SUFFIX="&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true&limit=0"

usage() {
    echo "Usage: $0 <accession>"
    echo "Example: $0 PRJEB47011"
    exit 1
}

# Check if the accession is provided
if [ -z "$1" ]; then
    usage
    exit
fi

URL="${URL_PREFIX}${1}${URL_SUFFIX}"

# Download URL to 1.tsv
curl -s "$URL" > $1.tsv