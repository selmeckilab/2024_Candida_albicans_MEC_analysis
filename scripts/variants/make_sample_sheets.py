#!/usr/bin/env python3
"""
Author : nancyscott <scot0854@umn.edu>
Date   : 2025-03-13
Purpose: Generate sample files for each patient series
"""

import argparse
from typing import NamedTuple
import pandas as pd

class Args(NamedTuple):
    """ Command-line arguments """
    file: str
    dir: str
    study: str


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Generate sample files for each patient series',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file',
                        metavar='FILE',
                        help='An excel or CSV file')
    
    parser.add_argument('dir',
                        metavar='DIR',
                        help='Output directory')
    
    parser.add_argument('study',
                        metavar='STUDY',
                        help='study to subset on')

    args = parser.parse_args()

    return Args(args.file, args.dir, args.study)


# --------------------------------------------------
def main() -> None:
    """ Read metadata, group series and print sample files """

    args = get_args()
    
    if args.dir.endswith('/'):
        dir = args.dir
    else:
        dir = args.dir + '/'
    
    metadata = pd.read_excel(args.file, index_col="sample")
    metadata = metadata[metadata["study"]==args.study]
    metadata = metadata[["series"]].dropna()
    
    to_write = metadata.groupby("series").groups
    
    for k, vals in to_write.items():
        filename = dir + k + "_samples.txt"
        with open(filename, 'w') as f:
            for v in vals:
                f.write(v + "\n")
    
    
# --------------------------------------------------
if __name__ == '__main__':
    main()
