# -*- coding: utf-8 -*-
"""
Created on Wed May  8 13:43:03 2024

@author: ivanp
"""

import os
import argparse
import logging
from typing import List, Tuple, Dict, Union

import numpy as np
import pandas as pd


from .const import (CONTIG_XMAP, CONTIG_QUERY,
                    CHROM_XMAP, CHROM_REFERENCE, CHROM_QUERY,
                    REF_TOL, CON_TOL, MOL_TOL, DIS_TOL)
from .funcs import calculate_telomere_lengths
from .utils import joinpaths

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)
logger = logging.getLogger("telompy")

__all__ = ["command_line_target", "validate_targets_target"]

parser = argparse.ArgumentParser(
    description="Extract lengths of telomeres from the de novo assembly of BNGO data")

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

# FUNCTIONAL ARGUMENTS
parser.add_argument("-a", "--arms", type=str, nargs="+", choices=['l', 'L', 'R', 'r'], default=["L", "R"],
                    help="Calculate telomeres on which arm of a chromosome - left (L) or right (R) or both (L R)")

parser.add_argument("-rt", "--ref_tol", type=int, default=REF_TOL,
                    help="Maximum number of unpaired labels on reference (AFTER LAST PAIR)")
parser.add_argument("-ct", "--con_tol", type=int, default=CON_TOL,
                    help="Maximum number of unpaired labels on contig (AFTER LAST PAIR)")
parser.add_argument("-mt", "--mol_tol", type=int, default=MOL_TOL,
                    help="Maximum number of unpaired labels on molecule (AFTER LAST PAIR)")
parser.add_argument("-dt", "--dis_tol", type=int, default=DIS_TOL,
                    help="Maximal distance from the first/last aligned label to the respective end of the chromosome")

# OUTPUT FOLDER
parser.add_argument("-o", "--output", default="telomere_lengths", type=str,
                    help="Output folder for extracted data")

# TODO - refactor this from constants
# ARGUMENTS FOR FUTURE-PROOFING
parser.add_argument("-cf", "--contig_xmap_format", type=str, default=CONTIG_XMAP,
                    help="Reconfigure contig path")
parser.add_argument("-mx", "--chrom_xmap", type=str, default=CHROM_XMAP,
                    help="Reconfigure master xmap")
parser.add_argument("-qc", "--contig_query_format", type=str, default=CONTIG_QUERY,
                    help="Reconfigure format of query cmap (molecule query)")
parser.add_argument("-mr", "--chrom_reference", type=str, default=CHROM_REFERENCE,
                    help="Reconfigure format of reference cmap (contig reference/chromosome query)")
parser.add_argument("-cq", "--chrom_query", type=str, default=CHROM_QUERY,
                    help="Reconfigure format of contig as a query")


def load_config(conf: pd.DataFrame) -> Union[None, pd.DataFrame]:
    if conf.shape[1] < 2:
        logger.error("Two columns must be in configuration file")
        return None
    conf[1] = conf.apply(
        lambda row: os.path.basename(row[0]) if pd.isna(row[1]) else row[1],
        axis=1)
    return conf


def target_from_config(args: Dict[str, List[str]]) -> Union[None, pd.DataFrame]:
    "constructs target from a given config file"
    if args.get("conf", None) is None:
        return None
    try:
        conf = pd.read_csv(args["conf"], header=None)
    except FileNotFoundError:
        logger.error("No configuration file found at '%s'", args["conf"])
        return None
    return load_config(conf)


def target_from_input(args: Dict[str, List[str]]) -> Union[None, pd.DataFrame]:
    "constructs target from -I -N options"
    if args["input"] is None:
        return None

    inputs = args["input"]
    names = args["name"]
    names = list() if names is None else names
    while len(names) < len(inputs):
        # ignoring this for type hinting
        # easier to check np.nan
        names.append(np.nan)  # type: ignore
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
            logger.error(
                "No file at %s - will exclude it from calculation", path)
    return new_targets


def validate_targets_target(input_args: Dict[str, List[str]]) -> List[Tuple[str, str]]:
    "validates target"

    target_conf = target_from_config(input_args)
    target_in = target_from_input(input_args)

    if target_conf is None and target_in is None:
        raise ValueError("Must provide input either via --input or via --conf")

    # create targets and validate them
    targets = target_conf if target_conf is not None else target_in

    # converted to tuple
    targets = [tuple(row) for row in targets.to_records(index=False)]  # pylint:disable=E1101

    targets = validate_targets(targets)

    return targets


def reconfigure_arms(arms=List[str]
                     ) -> List[str]:
    "reconfigures arums"
    return [{"L": "left", "R": "right"}[k] for k in list(set(x.upper() for x in arms))]


def command_line_target():
    "main function - target for CLI"
    args = vars(parser.parse_args())
    # arguments for input
    input_args = {k: v for k, v in args.items() if k in {
        "input", "name", "conf"}}

    # arguments for func
    func_args = {k: v for k, v in args.items() if k in
                 {"contig_xmap_format", "chrom_xmap",
                  "contig_query_format", "chrom_reference", "chrom_query"
                  "ref_tol", "con_tol", "mol_tol", "dis_tol",
                  "how"}}
    func_args["how"] = reconfigure_arms(args["arms"])

    targets = validate_targets_target(input_args)

    output_dir = args["output"]
    os.makedirs(output_dir, exist_ok=True)
    for path, name in targets:

        data = calculate_telomere_lengths(path, **func_args)
        output_path = joinpaths(output_dir, f"{name}.csv")

        # pd.concat(data, axis=0).to_csv(output_path, index=False)
        data.to_csv(output_path, index=False)
        logger.info("Saved telomere lengths at %s", output_path)
