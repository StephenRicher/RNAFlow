#!/usr/bin/env python3

""" Merge featureCount outputs to count matrix. """

import re
import sys
import logging
import argparse
import pandas as pd
from typing import List
from pathlib import Path


def main(infiles: List, samples: List) -> None:

    if len(samples) != len(infiles):
        # Raise warning if sample names are provided by not the correct number
        if samples:
            logging.error('Length of sample names must match length of infiles')
        # By default set sample name to filename (without path)
        samples = [str(Path(file)) for file in infiles]

    countMatrix = []
    for sample, file in zip(samples, infiles):
        countMatrix.append(pd.read_csv(
            file, names=['geneID', sample], header=0,
            index_col=0, usecols=[0,6], sep='\t', comment='#'))

    pd.concat(countMatrix, axis=1).to_csv(sys.stdout)


def parse_arguments():
    """ Parse command line arguments. """

    epilog='Stephen Richer, University of Bath, Bath, UK (sr467@bath.ac.uk)'
    parser = argparse.ArgumentParser(epilog=epilog, description=__doc__)
    parser.add_argument(
        'infiles', nargs='+',
        help='featureCounts results to merge')
    parser.add_argument(
        '--samples', nargs='+', default=[],
        help='Sample names for counts matrix - must be one for each infile')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    sys.exit(main(**vars(args)))
