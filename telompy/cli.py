# -*- coding: utf-8 -*-
"""
Created on Wed May  8 13:43:03 2024

@author: ivanp
"""

import os
import argparse
from typing import List, Tuple, Dict

import numpy as np
import pandas as pd

from .const import LOGGER as logger
from .const import CONTIG_PATH, QUERYCMAP_PATH, MASTER_XMAP, MASTER_REFERENCE
from .funcs import calculate_telomere_lengths
from .utils import joinpaths


__all__ = ["command_line_target", "validate_targets_target"]

parser = argparse.ArgumentParser(description="Extract lengths of telomeres from the de novo assembly of BNGO data")

# POINTERS TO FILES
parser.add_argument("-i", "--input", type=str, nargs="+", required=False,
                    help="folder(s) where BNGO de novo assembly output is stored")
parser.add_argument("-n", "--name", type=str, nargs="+", required=False,
                    help="prefix(es) of how the files will be stored")

parser.add_argument("-c", "--conf", required=False,
                    help="Path to the configuration file for a from-file extracting")
# MULTITHREADING
parser.add_argument("-t", "--threads", type=int, default=1,
                    help="Number of threads for parallel extraction")

# OUTPUT FOLDER
parser.add_argument("-o", "--output", default="telomere_lengths", type=str,
                    help="Output folder for extracted data")

# ARGUMENTS FOR FUTURE-PROOFING
parser.add_argument("-ct", "--contig_format", type=str, default=CONTIG_PATH,
                    help="Reconfigure contig path")
parser.add_argument("-mx", "--main_xmap", type=str, default=MASTER_XMAP,
                    help="Reconfigure master xmap")
parser.add_argument("-cq", "--querycmap_format", type=str, default=QUERYCMAP_PATH,
                    help="Reconfigure format of query cmap")

parser.add_argument("-mr", "--main_cmapr", type=str, default=MASTER_REFERENCE,
                    help="Reconfigure format of query cmap")


def target_from_input(args: dict):
    "constructs target from -I -N options"
    if args["input"] is None:
        return None

    inputs = args["input"]
    names = args["name"]
    names = list() if names is None else names
    while len(names) < len(inputs):
        names.append(None)
    return list(zip(inputs, names))


def target_from_config(args: dict):
    "constructs target from a given config file"
    if args["conf"] is None:
        return None
    try:
        conf = pd.read_csv(args["conf"], header=None)
    except FileNotFoundError:
        logger.error("No configuration file found at '%s'", args["conf"])
        return None

    if conf.shape[1] < 2:
        logger.error("Two columns must be in configuration file")
        return None
    # E1101 (no-member) : Instance of 'TextFileReader' has no 'replace'
    # this is a false positive
    conf = conf.replace({np.nan: None})  # pylint: disable=E1101
    return list(zip(conf[0], conf[1]))


def redefine_targets(targets: list) -> list:
    "where no name is provided, take the basename of folder"

    for i, (path, name) in enumerate(targets):
        if name is None:
            base_name = os.path.basename(path)
            targets[i] = (path, base_name)
    return targets


def validate_targets(targets: list) -> list:
    "validates targets"

    new_targets = list()
    for tuple_pair in targets:
        path, name = tuple_pair
        if os.path.isdir(path):
            logger.info("Found file at %s - will name it %s", path, name)
            new_targets.append((path, name))
        else:
            logger.error("No file at %s - will exclude it from calculation",path)
    return new_targets


def validate_targets_target() -> Tuple[List[Tuple[str, str]], str, Dict]:
    "validates target"
    args = vars(parser.parse_args())

    target_conf = target_from_config(args)
    target_in = target_from_input(args)
    print(args)
    if target_conf is None and target_in is None:
        raise ValueError("Must provide input either via --input or via --conf")

    # create targets and validate them
    targets = target_conf if target_conf is not None else target_in
    targets = validate_targets(targets)

    # pop output from arguments
    output_dir = args.pop("output")
    # pops the argument
    for _arg in ["conf", "input", "name"]:
        if _arg in args:
            del args[_arg]

    return targets, output_dir, args


def command_line_target():
    "main function - target for CLI"

    targets, output_dir, args = validate_targets_target()
    os.makedirs(output_dir, exist_ok=True)
    for path, name in targets:
        logger.info("Calculating telomere length for file found at %s", path)
        data = calculate_telomere_lengths(path, **args)
        output_path = joinpaths(output_dir, f"{name}.csv")
        # TODO - align with output of data
        pd.concat(data, axis=1).to_csv(output_path)
        logger.info("Saved telomere lengths at %s", output_path)


def command_line_targetO():
    "main function - target for CLI"
    args = vars(parser.parse_args())

    target_conf = target_from_config(args)
    target_in = target_from_input(args)

    output_dir = args.pop("output")

    # pops the argument
    for _arg in ["conf", "input", "name"]:
        if _arg in args:
            del args[_arg]

    if target_conf is None and target_in is None:
        raise ValueError("Must provide input either via --input or via --conf")

    os.makedirs(output_dir, exist_ok=True)
    # create targets
    targets = target_conf if target_conf is not None else target_in

    # redefine targets
    targets = redefine_targets(targets)

    # afterwards, validate them
    #targets = validate_targets(targets)

    # iterate through them - no need for multiprocessing
    for path, name in targets:
        logger.info("Calculating telomere length for file found at %s", path)
        data = calculate_telomere_lengths(path, **args)
        output_path = joinpaths(output_dir, f"{name}.csv")
        # TODO - align with output of data
        pd.concat(data, axis=1).to_csv(output_path)
        logger.info("Saved telomere lengths at %s", output_path)
