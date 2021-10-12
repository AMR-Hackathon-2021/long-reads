#!/usr/bin/env python
import sys
import pandas as pd
import argparse

# Program version
VERSION = '0.1'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Average coverage')
    parser.add_argument('--input', '-i', help="Input BED file",required=True )
    parser.add_argument('--output', '-o', help="Output table",required=False )
    parser.add_argument('--method', '-m', help="Method [default: mean]",default="mean" )
    parser.add_argument( '-s', '--separator', help="Output separator [default: tab]", default="\t" )
    parser.add_argument( '--header', help="Print header", action='store_true', default=False  )
    
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t', header=None, names=['chr', 'start', 'end', 'cov'])
    df['cov'] = df['cov'].astype(float)

    # Calculate average by chromosome
    try:
        df_avg = df.groupby('chr').agg({'cov': args.method})
    except ValueError:
        print("Error: method {} not supported".format(args.method))
        sys.exit(1)
        
    df_avg.reset_index(inplace=True)
    df_avg.columns = ['chr', 'cov']
    # Set column "chr" as index
    df_avg.set_index('chr', inplace=True)

    # Print output
    if args.output:
        df_avg.to_csv(args.output, sep=args.separator, header=args.header)
    else:
        df_avg.to_csv(sys.stdout, sep=args.separator, header=args.header)