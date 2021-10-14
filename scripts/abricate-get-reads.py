#!/usr/bin/env python
"""
Extract reads mapping to a ABRICATE AMR predictions.

Input: 
  + AMR predictions in tabular format
  + BAM file with reads mapped against the assembly

Output:
  + FASTQ file with reads mapping to the AMR predictions
"""

VERSION = "0.2"
AUTHOR = "Andrea Telatin"
LICENSE = "MIT"


import argparse
import sys
import os
import pandas as pd
import pysam

class Read(object):
    """
    A FASTQ read object, create as (name, seq, qual) or (name, comment, seq, qual)
    """

    def __init__(self, name, comment, seq, qual=None):
        self.name = name
        if qual:
            # read, comment, seq, qual
            if len(seq) != len(qual):
                raise ValueError("Sequence and quality length mismatch")
            self.comment = comment
            self.seq = seq
            self.qual = qual
        elif len(comment) == len(seq):
            # read, seq, qual
            self.comment = None
            self.seq = comment
            self.qual = seq
        else:
            raise ValueError("Invalid Read() from: {}, {}, {}".format(name, comment, seq))

    def __str__(self):
        comment = "\t" + self.comment if self.comment else ""
        return "@" + self.name + comment + "\n" + self.seq + "\n+\n" + self.qual

    # Length
    def __len__(self):
        return len(self.seq)
    
    # subscript sequence adn quality
    def __getitem__(self, index):
        return Read(self.name, self.seq[index], self.qual[index])
    
    def revcompl(self):
        """
        Reverse complement a sequence
        """
        complement = {
            'A': 'T', 
            'C': 'G', 
            'G': 'C', 
            'T': 'A'}
        return Read(self.name, "".join([complement[base] for base in self.seq[::-1]]), self.qual[::-1])

def eprint(*args, **kwargs):
    """
    Print to stderr
    """
    print(*args, file=sys.stderr, **kwargs)

def getBamAndTabFiles(directory):
    """
    Find the .bam and .tab files in a directory
    It's also grabbing the fasta files just in case...
    """
    files = [None, None, None]
    for f in os.listdir(directory):
        if f.endswith('.fasta'):
            if files[0] is None:
                files[0] = os.path.join(directory, f)
            else:
                # Being optional we just remove the FASTA from the list
                files[0] = None
        elif f.endswith('.bam'):
            if files[1] is None:
                files[1] = os.path.join(directory, f)
            else:
                raise Exception("Multiple .bam files found")
        elif f.endswith('.tab'):
            if files[2] is None:
                files[2] = os.path.join(directory, f)
            else:
                raise Exception("Multiple .tab files found")
            
    if files[1] is None or files[2] is None:
        raise Exception("ERROR: .bam or .tab files are missing in {}".format(directory))

    return files

def qualityStringFromQualities(qualities):
    """
    Convert a list of quality scores to a string
    """
    if qualities is None:
        raise ValueError("No quality scores")
    return "".join([chr(x + 33) for x in qualities])

def getReadsFromCoordinates(bamfile, chrname, start, end):
    """
    Extract reads from a BAM file given a region
    """
    
    bamiter = bamfile.fetch(chrname, start, end)
    c = 0

    for read in bamiter:
        c += 1

        # Skip not primary alignments
        if  read.is_secondary :
            continue

        name = read.query_name
        if start + read.infer_query_length() < end:
            # This read is not fully contained in the region
            continue

        comment = f"chr={read.reference_name};mapq={read.mapping_quality};"
        seq = read.get_forward_sequence()
        qual = read.get_forward_qualities()

        
        try:
 
            qual = qualityStringFromQualities(read.get_forward_qualities())
        except ValueError:
            raise ValueError("No quality scores: {}, {}".format(name, qual))
        except:
            raise Exception("Unexpected error: {}".format(sys.exc_info()[0]))
        
        try:
            fastqread = Read(name, comment, seq, qual)
            yield fastqread
        except ValueError:
            sys.stderr.write("WARNING: Read {} has different sequence and quality length\n".format(name))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    group_input = parser.add_argument_group('Input files')
    group_input.add_argument('-i', '--dir', required=False, help='Directory with assembly, bam file and Abricate file')
    group_input.add_argument('-a', '--amr', required=False, help='ABRICATE predictions in tabular format')
    group_input.add_argument('-b', '--bam', required=False, help='BAM file with reads mapped against the assembly')
    #group_input.add_argument('-g', '--genome', required=False, help='FASTA file with the assembly [optional]')

    group_output = parser.add_argument_group('Output')
    group_output.add_argument('-o', '--out-dir', required=False, help='Output directory (will split FASTQ files by gene)')
    group_output.add_argument('--format', choices=['gene', 'accession', 'product'], default="accession", help="Output filename [gene, product, accession]. Default: accession")
    group_output.add_argument('--min-len', dest="minlen",   help="Skip reads shorter than MINLEN [default: 500]", default=500)

    misc = parser.add_argument_group('Options')
    misc.add_argument( '--min-cov',  help='Minimum coverage in ABRICATE [90.0]', default=90.0, type=float)
    misc.add_argument( '--min-id',  help='Minimum identity quality in ABRICATE [94.0]', default=94.0, type=float)
    misc.add_argument( '--gene', help='Gene name to extract [None]', default=None)
    misc.add_argument('-v', '--verbose', action='store_true', help='Verbose mode')

    args = parser.parse_args()

    fastq_file, bam_file, amr_file = None, None, None

    # Check that --dir is set or --amr, --fastq and --bam are set
    if args.dir:
        fastq_file, bam_file, amr_file = getBamAndTabFiles(args.dir)
    elif args.amr and args.fastq and args.bam:
        amr_file = args.amr
        bam_file = args.bam
    else:
        sys.exit('ERROR: Please provide either --dir INPUTDIR or --amr ABRICATE and --bam BACKMAPPING')


    # Load abricate table
    try:
        amr_table = pd.read_csv(amr_file, sep='\t')
    except IOError:
        sys.exit('ERROR: Cannot open abricate table %s' % amr_file)

    # Check that required columns are in the amr table
    required_cols  = ["GENE", "ACCESSION", "PRODUCT", "%COVERAGE", "%IDENTITY"]
    for col in required_cols:
        if col not in amr_table.columns:
            sys.exit('ERROR: Column %s not in abricate table' % col)

    # Create outdir if not exists
    if args.out_dir:
        if not os.path.exists(args.out_dir):
            os.makedirs(args.out_dir)
    # Check BAM index 
    if not os.path.exists(bam_file + '.bai'):
        sys.exit('ERROR: Cannot open %s.bai, is the BAM file indexed?' % bam_file)

    # Load BAM file
    try:
        bamfile = pysam.AlignmentFile(bam_file, "rb")
    except IOError:
        sys.exit('ERROR: Cannot open %s' % bam_file)
    except ValueError:
        sys.exit('ERROR: Cannot open %s' % bam_file)


    # Remove predictions with coverage or identity below thresholds
    amr_table = amr_table[amr_table['%COVERAGE'] >= args.min_cov]
    amr_table = amr_table[amr_table['%IDENTITY'] >= args.min_id]

    # Remove predictions where PRODUCT does not contain args.gene
    if args.gene:
        amr_table = amr_table[amr_table['PRODUCT'].str.contains(args.gene)]

    # apply getReadsFromCoordinates to each row of the table
    for index, row in amr_table.iterrows():
        if args.verbose:
            eprint("Processing {}".format(row['PRODUCT']))
            
        if args.out_dir:
            outfile = os.path.join(args.out_dir, row[(args.format).upper()] + '.fastq')

            # open outfile for writing
            try:
                outfile = open(outfile, 'w')
            except IOError:
                sys.exit('ERROR: Cannot open %s' % outfile)
        else:
            # Set outfile = stdout
            outfile = sys.stdout

        chrname = row['SEQUENCE']
        start = row['START']
        end = row['END'] 
        reads = getReadsFromCoordinates(bamfile, chrname, start, end)
        genelen = end - start
        n = 0
        for r in reads:
            n += 1
            r.comment += 'len={};gene={};product={};accession={};coverage={};identity={};'.format(genelen, row['GENE'], row['PRODUCT'], row['ACCESSION'], row['%COVERAGE'], row['%IDENTITY'])
            if len(r) < int(args.minlen):
                continue
            print(r, file=outfile)
    