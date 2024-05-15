# -*- coding: utf-8 -*-
"""
Created on Wed May  8 10:14:03 2024

@author: ivanp
"""
import os
import time
import functools
import platform

import pandas as pd
from .const import LOGGER as logger

__all__ = ["func_timer", "joinpaths", "read_map_file"]


def func_timer(func):
    "timer for a function"
    def wrapper_func(*args, **kwargs):
        "plc"
        time_s = time.time()
        func_result = func(*args, **kwargs)
        time_e = time.time()
        delta = time_e - time_s

        # logger.info(f"'Operation done in {delta:.4f} seconds")
        logger.info("Operation done in %.2f seconds", delta)

        return func_result
    return wrapper_func

# FUNCTIONS TO NORMALIZE PATH FINDING

# TODO - add absolute path expander


def windows_normalizer(p):  # pylint:disable=C0103
    "FORWARDSLASH -> BACKSLASH"
    return p.replace("/", "\\")


def posix_normalizer(p):  # pylint:disable=C0103
    "BACKSLASH -> FORWARDSLASH"
    return p.replace("\\", "/")


@functools.lru_cache(maxsize=None)
def detect_normf():
    "detects platform for path correction"
    if platform.system() == "Windows":
        return windows_normalizer
    if platform.system() in {"Linux", "Darwin"}:
        return posix_normalizer
    raise OSError("Unknown operating system")


normfunc = detect_normf()


def joinpaths(*paths):
    "joins and normalizes paths"
    return os.path.join(*[normfunc(p) for p in paths])

# FUNCTIONS FOR READING XMAP PATHS


def read_map_file(path):
    """
    Reads a *MAP file. *MAP files (XMAP,CMAP,SMAP)
    are tab delimited files that start with lines prefixed with
    '#' (header lines) and contain one line prefixed with '#h'
    which gives out the names of the columns
    """

    with open(path) as _f:
        i = 0
        header = None
        for line in _f.readlines():
            if line.startswith("#"):
                i += 1
                if line.startswith("#h"):
                    header = [x.strip() for x in line.replace("#h ", "").split("\t")]
            else:
                break

    data = pd.read_csv(path, sep="\t", skiprows=i, header=None)
    try:
        assert data.shape[1] == len(header)
    except AssertionError:
        if header[0] == "#h":
            header.pop(0)
        assert data.shape[1] == len(header)  # pylint:disable=E1101
    data.columns = header
    return data
