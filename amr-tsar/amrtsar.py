"""
amrtsar takes uncorrected long reads and detected ARGS and associated plasmids

amrtsar takes uncorrected long reads and detected ARGS and associated plasmids.
"""
from util import check_reads_ok, run_amr_prediction, run_plasmid_prediction, make_plot
import logging
import time 
import meta
import argparse


epi = "Licence: " + meta.__licence__ +  " by " +meta.__author__ + " <" +meta.__author_email__ + ">"
logging.basicConfig()
log = logging.getLogger()

def main(args):
    check_reads_ok()
    run_amr_prediction()
    run_plasmid_prediction()
    make_plot()
    print('You got amr dude')

if __name__ == '__main__':
    start_time = time.time()
    log.setLevel(logging.INFO)
    desc = __doc__.split('\n\n')[0].strip()
    parser = argparse.ArgumentParser(description=desc,epilog=epi)
    
    parser.add_argument ('-v', '--verbose', action='store_true', default=False, help='verbose output')
    parser.add_argument('--version', action='version', version='%(prog)s ' + meta.__version__)
    parser.add_argument('r1', action='store', help='Reads (FASTQ)')
    args = parser.parse_args()
    main(args)
    if args.verbose: 
        log.setLevel(logging.DEBUG)
        log.debug( "Executing @ %s\n"  %time.asctime())    
    if args.verbose: 
        log.debug("Ended @ %s\n"  %time.asctime())
        log.debug('total time in minutes: %d\n' %((time.time() - start_time) / 60.0))    