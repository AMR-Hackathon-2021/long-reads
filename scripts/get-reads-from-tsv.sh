#!/bin/bash
# Download the E. colURL reads from EBI ENA
# Taking as input a TSV file from EBI ENA with the following columns:
#  0: study_accession	
#  1: sample_accession	
#  2: experiment_accession	
#  3: run_accession	
#  4: tax_id	
#  5: scientific_name	
#  6: fastq_ftp	
#  7: submitted_ftp	
#  8: sra_ftp
# The `get-url-list-from-ena` script can be used to download such TSV file

# Absolute full path to script dir
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR=$(dirname "$SCRIPT_DIR")

# PROGRAM DEFAULTS
# TSV file with the URLs
INPUT_FILE=$REPO_DIR/docs/datasets-PRJEB47011.tsv
OUTPUT_DIR=$REPO_DIR/datasets/ecoli-reads-PRJEB47011
VERBOSE=0
NUM_SAMPLES=0
QUERY_STRING=""
# Print help funciton
print_help() {
    echo "Download the E. coli reads from EBI ENA"
    echo "Usage: get-ecoli-reads.sh [-v] [-i <input_file>] [-o <output_dir>]"
    echo ""
    echo "  -v:     verbose"
    echo "  -q STR: Require STR in URL"
    echo "  -i TSV: input file with the URLs in EBI ENA tsv format [default: $INPUT_FILE]"
    echo "  -o DIR: output directory [default: $OUTPUT_DIR]"
    echo "  -n INT: Stop after downloading INT samples"
}

while getopts ":vhi:o:n:q:" opt; do
  case $opt in
    v)
      VERBOSE=1
      ;;
    i)
      INPUT_FILE=$OPTARG
      ;;
    o) 
      OUTPUT_DIR=$OPTARG
      ;;
    n)
      NUM_SAMPLES=$OPTARG
      ;;
    q)
      QUERY_STRING="$OPTARG"
      ;;
    h)
        print_help
        exit 0
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

if [[ ! -e "$INPUT_FILE" ]]; then
    echo "Input file $INPUT_FILE does not exist"
    exit 1
fi

# Create output dir
mkdir -p "$OUTPUT_DIR"

# Parse TSV file
TOTAL=0
while IFS=$'\t' read -r -a line; do
    # 0 study_accession	
    # 1 sample_accession	
    # 2 experiment_accession	
    # 3 run_accession	
    # 4 tax_id	
    # 5 scientific_name	
    # 6 fastq_ftp	
    # 7 submitted_ftp	
    # 8 sra_ftp

    
    # Skip header
    if [[ "${line[0]}" == "study_accession" ]]; then
        continue
    fi
    TOTAL=$((TOTAL+1))
    # Stop if NUM_SAMPLES is is > 0 and we have downloaded NUM_SAMPLES samples
    if [[ $NUM_SAMPLES -gt 0 && $TOTAL -gt $NUM_SAMPLES ]]; then
        echo Downloaded $TOTAL/$NUM_SAMPLES samples
        break
    fi
    
    # Get sample accession
    PROJECT_ID="${line[0]}"
 
    # Get sample name
    SAMPLE_ID="${line[1]}"

    # Get sample URL
    FTP="${line[7]}"
 
    echo "[$TOTAL] $SAMPLE_ID"
    # Split $FTP on the ";"
    IFS=';' read -ra ADDR <<< "$FTP"
    for URL in "${ADDR[@]}"; do
        # Skip if URL does not contain $QUERY_STRING
        if [[ $QUERY_STRING != "" && ! "$URL" =~ $QUERY_STRING ]]; then
            echo "Skipping $URL: $QUERY_STRING not matched"
            continue
        fi
        echo Downloading $URL
        # Download file from FTP url
        wget -P "$OUTPUT_DIR" "$URL"
    done
done < "$INPUT_FILE"

