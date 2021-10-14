#!/usr/bin/env python
"""
Andrea Telatin 2021 - AMR Hackathon 

Generate artificial taxonomy to label
sequences from a FASTA file.
"""
import sys
import os
import argparse
import gzip


def nodesHeader():
    """
    Return header for nodes.dmp file
    1|1|no rank||8|0|1|0|0|0|0|0||
    """
    return "\t|\t".join(["1", "1", "no rank", "", "8", "0", "1", "0", "0", "0", "0", "0", ""])  + "\t|"

def namesHeader():
    """
    Return header for names.dmp file
    """
    return "\t|\t".join(["1", "root", "", "scientific name"]) + "\t|"

def makeNode(tax_id, parent=1, rank="genus"):
    """
    0 tax_id					-- node id in GenBank taxonomy database
 	1 parent tax_id				-- parent node id in GenBank taxonomy database
 	2 rank					-- rank of this node (superkingdom, kingdom, ...) 
 	3 embl code				-- locus-name prefix; not unique
 	4 division id				-- see division.dmp file
 	5 inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
 	6 genetic code id				-- see gencode.dmp file
 	7 inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
 	8 mitochondrial genetic code id		-- see gencode.dmp file
 	9 inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
 	10 GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
 	11 hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
 	12 comments				-- free-text comments and citations
    """
    data = [
        tax_id, # 0 tax_id
        parent, # 1 parent tax_id
        rank,   # 2 rank
        "",    # 3 embl code
        "",   # 4 division id 
        "",  # 5 inherited div flag 
        "", # 6 genetic code id 
        "", # 7 inherited GC  flag 
        "", # 8 mitochondrial genetic code id 
        "", # 9 inherited MGC flag 
        "", # 10 GenBank hidden flag 
        "", # 11 hidden subtree root flag 
        ""] # 12 comments

    # convert every data element to string
    data = [str(x) for x in data]
    return "\t|\t".join(data) + "\t|"

def makeName(taxid, name):
    return "\t|\t".join([str(taxid), str(name), "", "scientific name"]) + "\t|"

def giHeader():
    return "\t".join(["accession", "accession.version", "taxid", "gi"]) 

def makeGi(taxid, name):
    return "\t".join([name, name, str(taxid), str(taxid)])

def eprint(*args, **kwargs):
    """
    Print to STDERR
    """
    print(*args, file=sys.stderr, **kwargs)
 


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
                        help='Sequences in FASTA format')
    parser.add_argument('-o', '--outdir', required=True, help='Output directory')
    parser.add_argument('-s', '--split', required=False, help='Split sequence name on this char')
    parser.add_argument('-f', '--field', required=False, help='Select this field on the name (requires --split)')
    parser.add_argument('-r', '--regex', required=False, help='n/a')
    
    parser.add_argument('--taxonomy-field', required=False, default="authority", help="Selector for names.dmp [default: authority]")
    args = parser.parse_args()

    if args.split and not args.field:
        parser.error("--field and --split are both required.")

    # Create output directory
    if not os.path.exists(args.outdir):
        try:
            os.makedirs(args.outdir)
        except OSError as exc: # Guard against
            eprint("Error creating output directory:", exc)
            sys.exit(1)
    
    # Create taxonomy file: names.dmp
    names_file = os.path.join(args.outdir, "names.dmp")
    nodes_file = os.path.join(args.outdir, "nodes.dmp")
    gi_file    = os.path.join(args.outdir, "nucl_gb.accession2taxid")
    seq_file =   os.path.join(args.outdir, "sequences.fa")

    # Create nodes and names file handles
    namesFh = open(names_file, "w")
    nodesFh = open(nodes_file, "w")
    giFh    = open(gi_file, "w")
    seqFh   = open(seq_file, "w")

    print( nodesHeader(), file=nodesFh)
    print( namesHeader(), file=namesFh)
    print( giHeader(), file=giFh)
    # Parse fasta
    seqId = 2

    for name, seq in fasta_iter(args.input):
        #>gb|JN967644|+|0-813|ARO:3002356|NDM-6 [Escherichia coli]
        # search for [Species name] in the name
        seqname = name.split(" ")[0].split("\t")[0]

        seqId += 1

        if args.split:
            seqname = seqname.split(args.split)[int(args.field) - 1]

        print(makeNode(seqId), file=nodesFh,)
        print(makeName(seqId, seqname), file=namesFh)
        print(makeGi(seqId, seqname), file=giFh)
        print(f">{seqname}|kraken:taxid|{seqId}\n{seq}", file=seqFh)
        
"""
The names.dmp file should have these columns with these delimiters:
"taxid\t|\tdisplay_name\t|\t-\t|\tscientific name\t|\n"
9606 | Homo sapiens |  | scientific name | 

The nodes.dmp file should have these columns with these delimiters:
"taxid\t|\tparent_taxid\t|\trank\t|\t-\t|\n"
1 | 1 | no ranx |
"""