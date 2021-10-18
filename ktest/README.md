### Exploratory analysis of kmer composition of reads and databases
#### using compressed sequences to reduce the noise of homopolymer errors.

* `cmers.nim` (compressed kmers) can generate a table of kmers present in a reference (e.g. AMR database) and
  then scan long reads for these kmers.


### Compiling

1. Install Nim (compiler and _nimble_ package manager) if not available in the system
1. Clone [SeqFu repository](https://github.com/telatin/seqfu2) for libraries to `$SEQFU_DIR`
2. From `$SEQFU_DIR` run `nimble install --depsOnly` to install the required dependencies
3. From this directory `nim c --threads:on -p:${SEQFU_DIR}/src  -p:${SEQFU_DIR}/src/lib cmers.nim`


### Usage

```
Compressed-mers

  A program to select long reads based on a compressed-mers dictionary

  Usage: 
  cmers scan [options] <DB> <FASTQ>...
  cmers make [options] <DB> 

  Make db options:
    -k, --kmer-size INT    K-mer size [default: 15]
    -o, --output-file STR  Output file [default: stdout]
    --subsample FLOAT  Keep only FLOAT% kmers [default: 100.0]
    --keep-multi INT       Keep kmers with multiple more than INT hits (when subsampling) [default: 0]
    --keep-single          Keep all kmers with single hit (when subsampling)

  Scanning options:
    -w, --window-size INT  Window size [default: 1500]
    -s, --step INT         Step size [default: 350]
    --min-len INT          Discard reads shorter than INT [default: 500]
    --min-hits INT         Minimum number of hits per windows [default: 50]
  
  Multithreading options:
    --pool-size INT        Number of sequences per thread pool [default: 1000]
    --max-threads INT      Maximum number of threads [default: 64]
    
    --verbose              Print verbose log
    --help                 Show help
```