# -*- coding: utf-8 -*-
"""

@author: ivanp
"""
import os

import logging
from typing import List, Set, Dict, Generator, Optional, Tuple

import numpy as np
import pandas as pd

from .utils import func_timer,  get_partial_alignment, get_last_labels, count_comment_lines

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)
logger = logging.getLogger("telompy_fandom")


# %% efficiently reading and processing molecules in BNX file


def bnx_generator(path: str, queried_molecules: Set[int]) -> Generator[List[str], None, None]:
    "yields molecules found in a given query in the form of a string list"
    molecule_collection = []  # refactor this

    molid = None
    with open(path, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue

            if line.startswith("0"):
                molid = line.split("\t")[1].strip()
                if molid not in queried_molecules:
                    molid = None
            if line.startswith("1") and molid:
                labels = line.strip().split("\t")[1:]
                labels.insert(0, molid)
                molecule_collection.append(labels)
                molid = None
                yield labels


@func_timer
def read_queried_molecules(path: str, queried_molecules: Set[str]) -> List[np.array]:
    """
    Given BNX path, and a set of molecule ids (as strings),
    returns a list of numpy array.

    Every numpy array is a representation of molecule where
    value at index 0 is the molecule ID, and value at the last index
    is the molecule length, so everything between is a position of label.

    One benefit to this kind of approach is that BNGO indexing starts at 1,
    and array indexing starts at 0, so to get a nucleotide position
    of a label we can simply use the label value, instead of label_value - 1
    """
    # 2.49 s ± 17.2 ms per loop (mean ± std. dev. of 7 runs, 1 loop each) 100000

    # 10000 random
    return [np.array(a, dtype=float) for a in
            bnx_generator(path=path, queried_molecules=queried_molecules)]
# %% efficiently reading and processing XMAP alignments


def fish_last_label(xmap_df: pd.DataFrame, last_label_mapper: Dict[int, int]) -> pd.DataFrame:
    """
    Extracts last label on alignment
    using knowledge about what last label for a chromosome is

    Then returns a dataframe with 4 columns:
        QryContigID - ID of molecule
        RefContigID - ID of chromosome
        Orientation - 1 for + and -1 for -
        QueryLabel - Label on molecule paired with the last label on the chromosome
    """
    lastlabel = xmap_df.copy()

    # extract last pair of labels
    # by design - refrence labels are ascending
    # so we extract last pair
    # - it corresponds to last aligned label on refrence
    lastlabel["LP"] = lastlabel.Alignment.str.split(
        ")(", regex=False).str[-1].str[:-1]
    lastlabel[["REF_LAB", "QueryLabel"]
              ] = lastlabel.LP.str.split(",", expand=True)

    lastlabel["RefContigID"] = lastlabel["RefContigID"].astype(int)
    lastlabel["REF_LAB"] = lastlabel["REF_LAB"].astype(int)
    lastlabel["QueryLabel"] = lastlabel["QueryLabel"].astype(int)

    lastlabel = lastlabel.loc[lastlabel.REF_LAB ==
                              lastlabel["RefContigID"].map(last_label_mapper)]
    lastlabel["Orientation"] = lastlabel["Orientation"].map({"+": 1, "-": -1})

    return lastlabel[["QryContigID", "RefContigID", "Orientation", "QueryLabel"]]


@func_timer
def read_fandom_xmap(xmap_path: str, last_label_mapper: Dict[int, int],
                     chunksize: int = 50000, ) -> pd.DataFrame:
    """
    Extracts last label on alignment
    using knowledge about what last label for a chromosome is

    Then returns a dataframe with 4 columns:
        QryContigID - ID of molecule
        RefContigID - ID of chromosome
        Orientation - 1 for + and -1 for -
        QueryLabel - Label on molecule paired with the last label on the chromosome

    For a file ~ 1.3GB :40.2 s ± 6.05 s per loop (mean ± std. dev. of 7 runs, 1 loop each)
    """
    skiprows = count_comment_lines(xmap_path)

    with pd.read_csv(xmap_path, chunksize=chunksize, sep="\t",
                     skiprows=skiprows, usecols=[1, 2, 7, 13],
                     names=["QryContigID", "RefContigID",
                            "Orientation", "Alignment"]
                     ) as reader:
        elements = [fish_last_label(chunk, last_label_mapper)
                    for chunk in reader]
    return pd.concat(elements, ignore_index=True)

# %% calculate telomere length using formula


@func_timer
def calculate_telomeres(telomeric_mols: List[np.array], telomeric_alns: np.array) -> np.array:
    """

    Calculates telomere length according to formula:

    If orientation is + (positive) then the formula is
        ML - PL
    If orientation is - then the formula is:
        PL
    where ML is the molecule length (last index) and PL is the position of
    a molecule label that is paired to the label of a given reference label
    found in telomeric_alns

    """
    result_array = list()
    # blind to sorting
    for molecule in telomeric_mols:
        # get alignment of specific molecule
        align = telomeric_alns[telomeric_alns[:, 0] == molecule[0]]

        # get length according to formula
        tlen = molecule[align[0, 3]]

        # if orientation is positive
        if align[0, 2] == 1:
            tlen = molecule[-1] - tlen

        # stores in a list
        result_array.append([align[0, 1], molecule[0], tlen])
    return np.array(result_array)


def calculate_telomere_lengths(xmap_path: str, bnx_path: str, cmap_path: str,
                               partial=True) -> pd.DataFrame:
    """
    Calculates telomere length. Requires an aligment file in XMAP format,
    an input file in BNX format, and a reference genome file in CMAP format.

    Returns a dataframe with 4 columns:
        RefContigID : corresponding to chromosome
        MoleculeID : corresponding to molecule ID
        TelomereLen: corresponding to length of telomere
        FullAlignment: 1 if from full alignment, 0 if from partial alignment

    FaNDOM aligns in two rounds. In first round, it aligns molecules completely,
    and then it aligns molecules partially.

    Molecules that are either unaligned OR
    have mean alignment score less than 5000/label OR
    total alignment length is lower than 25 kb OR
    alignment does not cover 80% of query length are targeted for partial alignments.

    These alignments are found in "_partial.xmap" right next to the file if
    --nopartial option is turned off. We can also calculate telomeres from the
    partial alignments. By default, we calculate telomeres from partial
    molecules, unless 'partial' is set to false, or there is no partial
    XMAP near the full XMAP

    Some notes:
        Molecules can also be in full alignments and in partial alignments.
        Preliminary testing showed that those molecules tend to have the same telomere
        length with low amount of discordance.

    """

    logger.info("Getting labels from reference")
    last_label_positions = get_last_labels(cmap_path)

    logger.info("Getting alignments from XMAP - Full alignments")
    telomeric_alns_full = read_fandom_xmap(xmap_path, last_label_positions)
    queried_molecules_f = set(telomeric_alns_full["QryContigID"].astype(str).unique())

    pxmap_path = get_partial_alignment(xmap_path)
    if partial and os.path.isfile(pxmap_path):
        logger.info("Getting alignments from XMAP - Partial alignments")
        telomeric_alns_part = read_fandom_xmap(pxmap_path, last_label_positions)
        queried_molecules_p = set(telomeric_alns_part["QryContigID"].astype(str).unique())
    else:
        logger.info("Skipping partial alignments")
        queried_molecules_p = set()
        telomeric_alns_part = pd.DataFrame()
    # join the two so I can extract molecules in one swell foop
    queried_molecules = queried_molecules_f.union(queried_molecules_p)

    logger.info("Getting molecules from BNX")
    telomeric_mols = read_queried_molecules(bnx_path, queried_molecules)

    # necessary to avoid index 0 out of bounds - happens when we have a molecule
    # that has no alignment
    telomeric_mols_full = [x for x in telomeric_mols if str(int(x[0])) in queried_molecules_f]
    telomeric_mols_part = [x for x in telomeric_mols if str(int(x[0])) in queried_molecules_p]

    logger.info("Calculating telomere lengths - Full alignments")
    telomeres_full = calculate_telomeres(telomeric_mols_full, telomeric_alns_full.values)
    telomeres_full = pd.DataFrame(data=telomeres_full, columns=[
        "RefContigID", "MoleculeID", "TelomereLen"]).sort_values("RefContigID")
    telomeres_full["FullAlignment"] = 1

    if partial and len(telomeric_alns_part) > 0:
        logger.info("Calculating telomere lengths - Partial alignments")
        telomeres_part = calculate_telomeres(telomeric_mols_part, telomeric_alns_part.values)
        telomeres_part = pd.DataFrame(data=telomeres_part, columns=[
            "RefContigID", "MoleculeID", "TelomereLen"]).sort_values("RefContigID")
        telomeres_part["FullAlignment"] = 0
    else:
        telomeres_part = pd.DataFrame(columns=telomeres_full.columns)

    telomeres = pd.concat([telomeres_full, telomeres_part], ignore_index=True)
    return telomeres.astype(int)

# %% ALIGNMENT STATISTICS - NOT PRESENT IN BIONANO!!!


def get_rough_xmap_stats(xmap_path: str, suffix: str = "") -> Tuple[Dict[str, float], Set[int]]:
    """
    Gets rough statistics and returns unique molecules that align to the reference
    Statistics include:
        Length of queried olecule (mean and std)
        Confidence of alignnment (mean and std)
        Number of molecules
        Number of unique molecules
    """
    if xmap_path is None:
        return dict(), set()

    skiprows = count_comment_lines(xmap_path)

    with pd.read_csv(xmap_path, chunksize=50000, sep="\t", skiprows=skiprows,
                     usecols=[1, 8, 10], names=["QryContigID", "MoleculeConfidence", "QryLen"]) as reader:
        elements = [chunk.astype(int) for chunk in reader]

    full = pd.concat(elements, ignore_index=True)

    stats_dict = dict()

    stats_dict[f"qrylen_mean_{suffix}"] = full["QryLen"].mean()
    stats_dict[f"qrylen_std_{suffix}"] = full["QryLen"].std()

    stats_dict[f"conf_mean_{suffix}"] = full["MoleculeConfidence"].mean()
    stats_dict[f"conf_std_{suffix}"] = full["MoleculeConfidence"].std()

    unique_mols = set(full["QryContigID"])

    stats_dict[f"molecules_{suffix}"] = len(full)
    stats_dict[f"unique_mols_{suffix}"] = len(unique_mols) / len(full)
    return stats_dict, unique_mols


def get_xmap_statistics(xmap_path: str, pxmap_path: Optional[str] = None) -> Dict[str, float]:
    """
    Gets rough statistics for two rounds of alignment via FaNDOM.
    Statistics include:
        Length of queried olecule (mean and std)
        Confidence of alignnment (mean and std)
        Number of molecules
        Number of unique molecules
    Second round is partial alignment, and additional statistics include:
        How many molecules are in both alignment ('joint_mols')
        Jaccard similarity of molecules of both alignments ('jaccard_btw_alns')
    """
    if pxmap_path is None or not os.path.isfile(pxmap_path):
        logger.error("Partial alignment not provided...inferring from path")
        pxmap_path = get_partial_alignment(xmap_path)
        if not os.path.isfile(pxmap_path):
            logger.error("Partial alignment not found...continuing without it")
            pxmap_path = None
        else:
            logger.info("Partial alignment inferred from path")

    name = os.path.basename(xmap_path).replace(".xmap", "")
    stats_full, mols_full = get_rough_xmap_stats(xmap_path, "full_alignment")
    stats_part, mols_part = get_rough_xmap_stats(pxmap_path, "partial_alignmet")

    stats = stats_full.copy()
    stats.update(stats_part)
    if mols_part:
        bothmols = len(mols_full.intersection(mols_part))
        jacc = len(mols_full.intersection(mols_part))/len(mols_full.union(mols_part))
        stats["joint_mols"] = bothmols
        stats["jaccard_btw_alns"] = jacc
    for k, v in stats.items():
        logger.critical("%s - STATS: '%s': %f", name, k, v)

    return stats
