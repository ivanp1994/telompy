# -*- coding: utf-8 -*-
"""
Created on Wed May  8 13:43:03 2024

@author: ivanp
"""
# TODO - test branch on the Padobran
import os
import argparse
import logging
from typing import List, Tuple, Dict, Union

import numpy as np
import pandas as pd


from .const import CONTIG_PATH, QUERYCMAP_PATH, MASTER_XMAP, MASTER_REFERENCE
from .funcs import calculate_telomere_lengths
from .utils import joinpaths

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)
logger = logging.getLogger("telompy")

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

parser.add_argument("-gs", "--gap_size", type=int, default=100000,
                    help="Gap size for calculation of Telomere Offset")

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


def load_config(conf: pd.DataFrame) -> Union[None, pd.DataFrame]:
    if conf.shape[1] < 2:
        logger.error("Two columns must be in configuration file")
        return None
    conf[1] = conf.apply(
        lambda row: os.path.basename(row[0]) if pd.isna(row[1]) else row[3],
        axis=1)


def target_from_config(args: Dict[str, str]) -> Union[None, pd.DataFrame]:
    "constructs target from a given config file"
    if args.get("conf", None) is None:
        return None
    try:
        conf = pd.read_csv(args["conf"], header=None)
    except FileNotFoundError:
        logger.error("No configuration file found at '%s'", args["conf"])
        return None
    return load_config(conf)


def target_from_input(args: Dict[str:str]) -> Union[None, pd.DataFrame]:
    "constructs target from -I -N options"
    if args["input"] is None:
        return None

    inputs = args["input"]
    names = args["name"]
    names = list() if names is None else names
    while len(names) < len(inputs):
        names.append(np.nan)
    return load_config(pd.DataFrame([inputs, names]).T)


def validate_targets(targets: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    "validates targets"

    new_targets = list()
    for tuple_pair in targets:
        path, name = tuple_pair
        if os.path.isdir(path):
            logger.info("Found file at %s - will name it %s", path, name)
            new_targets.append((path, name))
        else:
            logger.error("No file at %s - will exclude it from calculation", path)
    return new_targets


def validate_targets_target(input_args: Dict[str, str]) -> List[Tuple[str, str]]:
    "validates target"

    target_conf = target_from_config(input_args)
    target_in = target_from_input(input_args)

    if target_conf is None and target_in is None:
        raise ValueError("Must provide input either via --input or via --conf")

    # create targets and validate them
    targets = target_conf if target_conf is not None else target_in
    targets = validate_targets(targets)

    return targets


def command_line_target():
    "main function - target for CLI"
    args = vars(parser.parse_args())
    # arguments for input
    input_args = {k: v for k, v in args.items() if k in {"input", "name", "conf"}}

    # arguments for func
    func_args = {k: v for k, v in args.items() if k in
                 {"gap_size", "contig_format", "main_xmap", "querycmap_format", "main_cmapr"}}

    targets = validate_targets_target(input_args)

    output_dir = args["output"]
    os.makedirs(output_dir, exist_ok=True)
    for path, name in targets:
        logger.info("Calculating telomere length for file found at %s", path)
        data = calculate_telomere_lengths(path, **func_args)
        output_path = joinpaths(output_dir, f"{name}.csv")

        pd.concat(data, axis=0).to_csv(output_path)
        logger.info("Saved telomere lengths at %s", output_path)
