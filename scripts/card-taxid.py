#!/usr/bin/env python
"""
Andrea Telatin 2021 - AMR Hackathon 

Add taxid to the CARD AMR sequences.

Input:
    - CARD sequences in FASTA format
    - NCBI taxonomy (names.dmp)

Output:
    - CARD sequences with taxid added to the header
"""
import sys
import os
import argparse
import gzip
import re
def eprint(*args, **kwargs):
    """
    Print to STDERR
    """
    print(*args, file=sys.stderr, **kwargs)

def loadTaxonomy(taxonomy_file, label="authority") -> dict:
    """
    Load NCBI taxonomy and return a "Species name" -> "Taxid" dictionary.
    NOTE: This will use a species name the first two words. Collisions are expected
    and simply overwritten.
    """
    taxonomy = {}
    with open(taxonomy_file, 'r') as f:
        for line in f:

            #562	|	Escherichia coli (Migula 1895) Castellani and Chalmers 1919	|		|	authority	|
            fields = line.split('|')

            # Remove leading and trailing spaces from all fields
            fields = [field.strip() for field in fields]
            if fields[3] == label:
                # split species after the second space
                species = " ".join(fields[1].split(" ")[:2])
                taxonomy[species] = fields[0]
    return taxonomy


# FASTA File parser generator supporting gzipped files
def fasta_iter(path):
    """
    Given a fasta file, yield tuples of header, sequence
    """
    seqName = None
    with (gzip.open if path.endswith('.gz') else open)(path, 'rt') as fasta:
        for line in fasta:
            if line.startswith('>'):
                if seqName is not None:
                    yield seqName, prefix
                seqName = line[1:].rstrip()
                prefix = ""
            else:
                prefix += line.rstrip()
    yield seqName, prefix
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', required=True,
                        help='CARD sequences in FASTA format')
    parser.add_argument('-t', '--taxonomy', required=False,
                        help='NCBI taxonomy (names.dmp)')
    parser.add_argument('-o', '--output', required=False,
                        help='CARD sequences with taxid added to the header')
    parser.add_argument( '--taxonomy-field', required=False, default="authority", help="Selector for names.dmp [default: authority]")
    parser.add_argument('-f', '--fallback', dest="defTaxId", help="Fallback taxid if not found in taxonomy", default="32630")
    args = parser.parse_args()


    # Load taxonomy
    taxonomy = loadTaxonomy(args.taxonomy)    
    eprint("Loaded taxonomy:", len(taxonomy))

  
    # Parse fasta
    for name, seq in fasta_iter(args.input):
        #>gb|JN967644|+|0-813|ARO:3002356|NDM-6 [Escherichia coli]
        # search for [SPecies Name] in the name
        seqname = name.split(" ")[0].split("\t")[0]
        regex = r'\[([^\]]+)\]'
        match = re.search(regex, name)
        species = ""
        comment = ""
        taxid = args.defTaxId
        if match:
            species = match.group(1)
            if species in taxonomy:
                taxid = taxonomy[species] 

        if len(species) > 0:
            comment = " " + species
        print(">" + seqname + "|kraken:taxid|" + taxid + comment)
        print(seq)