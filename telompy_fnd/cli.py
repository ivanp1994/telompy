# -*- coding: utf-8 -*-
"""
Created on Fri May 24 18:38:50 2024

@author: Ivan
"""


import os
import argparse
import logging
import multiprocessing
from typing import List, Tuple, Dict, Union, Optional

import numpy as np
import pandas as pd

from .funcs import calculate_telomere_lengths, get_xmap_statistics
from .utils import joinpaths, normfunc


logger = logging.getLogger("telompy_fandom")

parser = argparse.ArgumentParser(description="Extract lengths of telomeres from the alignment of FaNDOM app")

# POINTERS TO FILES
parser.add_argument("-bnx", "--bnx", type=str, nargs="+", required=False,
                    help="BNX file (query for alignment)")
parser.add_argument("-xmap", "--xmap", type=str, nargs="+", required=False,
                    help="XMAP file (the output of FaNDOM alignment)")
parser.add_argument("-ref", "--ref", type=str, nargs="+", required=False,
                    help="Reference (CMAP) file")

parser.add_argument("-c", "--conf", required=False,
                    help="Path to the configuration file for a from-file extracting")

# OUTPUTS
parser.add_argument("-n", "--name", type=str, nargs="+", required=False,
                    help="prefix(es) of how the files will be stored")

parser.add_argument("-o", "--output", default="telomere_lengths",
                    help="Output folder for extracted data")

# MULTITHREADING
parser.add_argument("-t", "--threads", type=int, default=1,
                    help="Number of threads for parallel extraction")


# %% funcs

def load_config(conf: pd.DataFrame, ref: Optional[str] = None) -> Union[None, pd.DataFrame]:
    """
    loads configuration file
    """
    if conf.shape[1] < 2:
        logger.error("At least two columns must be in the configuration file ")
        return None
    # if this is the case then we have no reference
    # we take it from args['ref']
    if conf.shape[1] == 2:
        if ref is None:
            logger.error("No reference files provided")
            return None
        if len(ref) > 1:
            logger.error("If reference is passed with --conf, then only one file is permitted")
            return None
        conf[2] = ref[0]

    # if this is true we have no names column
    # we take it from 0 th column
    if conf.shape[1] == 3:
        conf[3] = conf[0].apply(lambda x: os.path.basename(x).replace(".xmap", ""))

    # we now check for nans
    if conf.shape[1] == 4:
        # xmap bnx ref names
        # do check for names

        conf[3] = conf.apply(
            lambda row: os.path.basename(row[0]).replace('.xmap', '') if pd.isna(row[3]) else row[3],
            axis=1
        )

    return conf


def target_from_config(args: Dict[str, str]) -> Union[None, pd.DataFrame]:
    "constructs target from a given config file"
    if args.get("conf", None) is None:
        return None
    try:
        conf = pd.read_csv(args["conf"], header=None)
    except FileNotFoundError:
        logger.error("No configuration file found at '%s'", args["conf"])
        return None
    return load_config(conf, args.get("ref", None))


def target_from_input(args: Dict[str, str]) -> Union[None, pd.DataFrame]:
    "constructs target from -bnx -xmap -ref options"

    if args.get("bnx", None) is None or args.get("xmap", None) is None:
        return None

    bnx_files = args.get("bnx", None)
    xmap_files = args.get("xmap", None)
    names = args.get("name", None)
    ref = args.get("ref", None)

    # process names
    names = list() if names is None else names
    while len(names) < len(bnx_files):
        names.append(np.nan)

    # check parity between molecules and aligments
    if len(bnx_files) != len(xmap_files):
        raise ValueError("Uneven numbers of BNX and XMAP files")

    # check parity between references and bnx files
    if ref is None:
        raise ValueError("No reference provided")

    if len(ref) == 1:
        ref = ref*len(bnx_files)
    else:
        if len(ref) < len(xmap_files):
            raise ValueError("Uneven references - pass either 1 or N references ")

    return load_config(pd.DataFrame([xmap_files, bnx_files, ref, names]).T, None)

# validate targets


def validate_targets(targets: List[Tuple[str, str, str, str]]) -> List[Tuple[str, str, str, str]]:
    "checks if the particular file exists and avoids it if it doesnt"
    new_targets = list()
    for target in targets:
        if not os.path.isfile(target[0]):
            logger.error("No XMAP file at %s - will exclude it from calculation", target[0])
            continue
        if not os.path.isfile(target[1]):
            logger.error("No BNX file at %s - will exclude it from calculation", target[1])
            continue
        if not os.path.isfile(target[2]):
            logger.error("No CMAP file at %s - will exclude it from calculation", target[2])
            continue
        new_targets.append(target)
    if not new_targets:
        logger.error("No valid targets")
        raise ValueError("No valid files, check error log")
    return new_targets


def validate_targets_target(args: Dict[str, str]) -> List[Tuple[str, str, str, str]]:
    "validates target files"

    target_conf = target_from_config(args)
    target_in = target_from_input(args)

    if target_conf is None and target_in is None:
        raise ValueError("Must provide input either via --input or via --conf")

    # create targets and validate them
    targets = target_conf if target_conf is not None else target_in

    # this is a pandas dataframe - we need to validate it

    output_dir = args.get("output", None)
    os.makedirs(output_dir, exist_ok=True)
    targets[3] = targets[3].apply(lambda x: joinpaths(output_dir, x))
    # normalize pathing for the entire dataset
    for _c in targets.columns:
        targets[_c] = targets[_c].apply(normfunc)

    # turn a target into a list of tuples
    targets = [tuple(row) for row in targets.to_records(index=False)] # pylint:disable=E1101

    # check if they're real
    targets = validate_targets(targets)

    return targets

# functions compatible with multiprocessing


def calculate_telomere_lengths_unpacker(xmap_path: str, bnx_path: str, cmap_path: str, output_file: str) -> None:
    "unpacks the calculation of telomeres - saving to file and returning stats"

    if not output_file.endswith(".csv"):
        output_file = output_file + ".csv"

    out_df = calculate_telomere_lengths(xmap_path, bnx_path, cmap_path)
    out_df.to_csv(output_file, index=False)
    get_xmap_statistics(xmap_path)


def calculate_telomere_lengths_wrapper(args: Tuple[str, str, str, str]) -> None:
    "wrapper around function - for multiprocessing"
    return calculate_telomere_lengths_unpacker(*args)

# cli target


def command_line_target() -> None:
    """
    target for command line arguments
    """
    args = vars(parser.parse_args())

    threads = args.pop("threads")
    targets = validate_targets_target(args)

    if threads == 1:
        for target in targets:
            calculate_telomere_lengths_wrapper(target)
    else:
        with multiprocessing.Pool(threads) as pool:
            # Use map to apply the function to each element in the targets list
            pool.map(calculate_telomere_lengths_wrapper, targets)