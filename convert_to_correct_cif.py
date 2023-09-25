#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 11:44:21 2023

@author: mojtaba
"""

import numpy as np
from ccdc import io
import argparse
from multiprocessing import Pool


def conver_to_csd_cif(input_file):
    crystal_reader = io.CrystalReader(input_file)
    with io.CrystalWriter(input_file, append=False) as crystal_writer:
            crystal_writer.write(crystal_reader[0])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "structures",
        type=str,
        nargs="+",
        help="CIF or RES files containing input structures."
             )
    parser.add_argument(
    "-j",
    "--num_threads",
    type=int,
    help="Number of parallel threads to use when doing the comparisons.",
    default=1
    )



    args = parser.parse_args()
    structures = [item for item in args.structures]
    my_func = np.vectorize(conver_to_csd_cif)
    # map(conver_to_csd_cif, structures)
    my_func(structures)

    # with Pool(int(args.num_threads)) as pool:
        # pool.map(conver_to_csd_cif, structures)

if __name__ == "__main__":
    main()
