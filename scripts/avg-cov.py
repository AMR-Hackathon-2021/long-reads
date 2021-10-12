#!/usr/bin/env python
import os, sys
import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Average coverage')
    parser.add_argument('--input', '-i', required=True, help='Input file')
    parser.add_argument('--output', '-o', required=False, help='Output file')
    parser.add_argument('--method', '-m', default="mean", help='Method')
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t', header=None, names=['chr', 'start', 'end', 'cov'])
    df['cov'] = df['cov'].astype(float)

    # Calculate average by chromosome
    df_avg = df.groupby('chr').agg({'cov': args.method})
    df_avg.reset_index(inplace=True)
    df_avg.columns = ['chr', 'cov']
    print(df_avg)