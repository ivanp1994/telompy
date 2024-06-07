# -*- coding: utf-8 -*-
"""
Created on Wed May  8 10:11:51 2024

@author: ivanp
"""
import os
import json
import logging
from typing import Iterable, Callable, Union, Literal, Optional, List
from functools import partial
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import stats

from .utils import read_map_file, joinpaths, func_timer
from .const import (CONTIG_XMAP, CONTIG_QUERY,
                    CHROM_XMAP, CHROM_REFERENCE, CHROM_QUERY,
                    REF_TOL, CON_TOL, MOL_TOL, DIS_TOL,
                    INFREP)

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)
logger = logging.getLogger("telompy")
pd.options.mode.chained_assignment = None  # default='warn'

# %% FUNCTIONS FOR FISHING FIRST/LAST LABELS ON MASTER ASSEMBLY


def fish_starting_aligned_label(xmap_df: pd.DataFrame,
                                how: Literal["left", "right"]) -> pd.DataFrame:
    """
    From a given xmap_df we find the either the
    last or the first aligned reference label
    and then we filter out contigs ('QryContigID')
    that do not align to said last aligned reference.

    If there are N chromosomes (referent contigs - 'RefContigID'), there
    should be M (M>=N) rows in the resulting output

    The first aligned reference label corresponds to the 'left' arm of
    chromosome (the 'left' telomere) and the last aligned reference label
    corresponds to the 'right' arm of chromosome (the 'right' telomere)
    """
    aligned = xmap_df.copy()

    # depending on the telomere arm we want either
    # first aligned label ('left') OR
    # last aligned label ('right')

    if how == "left":
        idx = 0
    if how == "right":
        idx = -1
    aligned["SiteID"] = aligned["Alignment"].str[1:-
                                                 1].str.split(")(", regex=False).str[idx]

    # get the first element of the pair - corresponding to the label on the reference
    aligned["SiteID"] = aligned["SiteID"].str.split(",").str[0].astype(int)

    # now grouping by chromosomes find the first (smallest) label
    # or the last (largest) label
    if how == "left":
        aligned = aligned.groupby("RefContigID").apply(
            lambda one_chrom: one_chrom[one_chrom["SiteID"] == one_chrom["SiteID"].min()])
    if how == "right":
        aligned = aligned.groupby("RefContigID").apply(
            lambda one_chrom: one_chrom[one_chrom["SiteID"] == one_chrom["SiteID"].max()])
    return aligned.reset_index(drop=True)


def fish_offsets(aligned: pd.DataFrame, reference_cmap: pd.DataFrame, how: Literal["left", "right"]):
    """
    Given an XMAP file in aligned that only has query contigs aligned
    to either the first ('left') or last ('right')
    aligned reference label and a CMAP file in reference_cmap this
    function finds the difference between the first/last label on reference and
    the first/last aligned label on the reference in the number of labels.

    """
    aligned_positions = aligned[["RefContigID",
                                 "SiteID", "QryContigID"]].drop_duplicates()
    aligned_positions = pd.merge(left=aligned_positions,
                                 right=reference_cmap[[
                                     "CMapId", "SiteID", "Position"]],
                                 left_on=["RefContigID", "SiteID"],
                                 right_on=["CMapId", "SiteID"],
                                 )[["RefContigID", "QryContigID", "SiteID", "Position"]]

    if how == "left":
        labels = reference_cmap[reference_cmap.SiteID ==
                                1][["CMapId", "SiteID", "Position"]]
    if how == "right":
        labels = reference_cmap[reference_cmap.SiteID == reference_cmap.NumSites][[
            "CMapId", "SiteID", "Position"]]

    labels.columns = ["RefContigID", "FirstID",
                      "BoundReferencePosition"]  # <-here

    # merge the two
    aligned_positions = pd.merge(left=aligned_positions,
                                 right=labels,
                                 left_on="RefContigID",
                                 right_on="RefContigID",
                                 how="left")

    aligned_positions["OffsetLabel"] = abs(
        aligned_positions["SiteID"] - aligned_positions["FirstID"])
    return aligned_positions


def fish_paired_label(path: str, how: Literal["left", "right"],
                      chrom_xmap: str = CHROM_XMAP,
                      chrom_reference: str = CHROM_REFERENCE,) -> pd.DataFrame:
    """
    Given a path to the BioNano denovo assembly, this function finds XMAP
    for assembly of contigs-to-chromosome, then finds alignments that
    align to either the last (the largest)
    or the first (the smallest) label on the reference.

    The last label corresponds to the 'right' arm of the chromosome (the 'right' telomere)
    and the first label corresponds to the 'left' arm of the chromosome (the 'left' telomere)

    In a situation where the last aligned label on the reference is not the
    last label on reference, an additional column called 'OffsetLabel' is found.
    This contains the difference between last/first label and last/firs aligned label
    in the number of labels between the two.
    """
    master_xmap = read_map_file(joinpaths(path, chrom_xmap))
    master_cmap = read_map_file(joinpaths(path, chrom_reference))

    # find either FIRST or LAST reference label on the master xmap
    # filtering out only contigs that align to it

    aligned_pair = fish_starting_aligned_label(master_xmap, how)

    # find the distance between the last aligned reference label
    # and the last reference label
    # if a contig maps onto the last label then the OffSet is 0
    aligned_offsets = fish_offsets(aligned_pair, master_cmap, how)

    # merge them
    aligned_pair = pd.merge(left=aligned_pair,
                            right=aligned_offsets[[
                                "QryContigID", "RefContigID", "OffsetLabel", "BoundReferencePosition"]],
                            left_on=["QryContigID", "RefContigID"],
                            right_on=["QryContigID", "RefContigID"])
    aligned_pair["TelomereArm"] = how

    return aligned_pair

# %% FUNCTIONS FOR CALCULATING TELOMERE LENGTHS


def get_first_query_bound_reference(align_str: str, how: Literal["left", "right"]) -> int:
    """
    An XMAP dataframe contains column Alignment which is a list of tuples of the type
    >>>> (ref_label_A,query_label_A)(ref_label_B,query_label_B)
    Which is pair of aligned labels.

    By convention, the first elements of a pair are in a ascending order
    (ref_label_B>ref_label_A, etc). Depending on the arm of chromosome we want,
    we either find the first ('left') pair or the last ('right') pair.
    We return the second element of the pair which corresponds to label on
    the query ('contig') that is paired with the bound label

    Example:
    align_str = '(33,17781)(34,17780)(35,17779)(36,17778)(37,17777)(38,17776)'

    If we want left (the first element) we take the pair (33,17781) and return
    17781.
    If we want right (the last element) we take the pair (38,17776) and return
    17776.

    Note that the labels in query (the second pair) can be either descending
    or ascending depending on the orientation
    """
    if how == "left":
        idx = 0
    if how == "right":
        idx = -1

    align_list = [x.split(",") for x in align_str[1:-1].split(")(")]
    return int(align_list[idx][1])


def get_number_of_unpaired_contig_labels(path: Union[str, pd.DataFrame], alignment: str, qrycontigid: int,
                                         orientation: Literal["-", "+"], telarm: Literal["left", "right"],
                                         master_orientation: Literal["-", "+"],
                                         chrom_query: str = CHROM_QUERY) -> int:
    """
    This function returns number of unpaired labels on the contig (query).
    The contig is a query in a contig-to-reference assembly.
    XMAP does not contain the information about the total number of labels on the
    query, so we need to open up master_query of that contig.


    The number of unpaired labels on query depends on whether we want
    AFTER (telarm 'right') or BEFORE (telarm 'left') the first pair on the reference.

    IMPORTANT NOTES:
    (1) The number of unpaired labels is calculated AFTER/BEFORE the LAST/FIRST paired
        label on the reference - so if there are unpaired labels on the query
        between the first and last pair, we don't count them.
    For example:
        C: ------|------|-------|--------
        M: -|---|------|---|----|-------

        In the above case we have one label in M before first pair, but also
        one unpaired label after first label. We only count one.
    (2) The number of unpaired labels on the query with respect to paired label
        on reference. We ignore unpaired labels of contig to molecule.
    For example:
        R: [-------------|--------------|
        C: -----|---|----|---------------

        In the above case, we have two labels on contig unpaired to reference.
        We ignore those two.

    DETAILS:

        # when we have left, we have the following scenarios
        # R : [------------|-------------
        #                  1
        # C : [---|---|----|-----|----|--
        # +       1   2    3     4    5
        # -       5   4    3     2    1

        # In the case of positive orientation, labels will be equal to LABELPOS - 1
        # In the case of negative orientation, labels will be equal to N - LABELPOS
        # where LABELPOS is the label that aligns to reference
        # and N is equal to total number of labels

        # when we have right, we have the following scenarios
        # R : --------------|--------------]
        #                   N
        # C : ----|----|----|-----|-----|--]
        # +       1    2    3     4     5
        # -       5    4    3     2     1

        # In the case of positive orientation, labels will be equal to N-LABELPOS
        # In the case of negative orientation, labels will be equal to LABELPOS-1

        # to sum up
        # LABELPOS-1 | (telarm=left,orientation=+) | (telarm=right,orientation=-)
        # N-LABELPOS | (telarm=left,orientation=-) | (telarm=right,orientation=+)


        # if we map {left:1,right:-1} , {+:1,-:-1}
        # LABELPOS-1 | (telarm=1,orientation=1) | (telarm=-1,orientation=-1)
        # N-LABELPOS | (telarm=1,orientation=-1) | (telarm=-1,orientation=1)

        # if we get mult
        # LABELPOS-1 | telarm x orientation = 1
        # N-LABELPOS | telarm x orientation -1
    """
    # in this case - telarm does not do anything
    # last_label_query = int(alignment[1:-1].split(",")[1])
    # buuut it does when I reuse this for the assembled contig -.-
    last_label_query = get_first_query_bound_reference(alignment, how=telarm)

    map_orient = {"+": 1, "-": -1}[orientation]
    map_telarm = {"left": 1, "right": -1}[telarm]
    map_morient = {"+": 1, "-": -1}[master_orientation]

    mult = map_orient * map_telarm * map_morient

    # if mult is positive we return last label
    if mult == 1:
        return last_label_query-1

    # if mult is negative, we must find N where N is the number of labels
    if isinstance(path, str):
        _cmap_contig_ref = read_molecules_cmap(
            joinpaths(path, chrom_query), [qrycontigid])

    # if the path is an already loaded CMAP
    elif isinstance(path, pd.DataFrame):
        _cmap_contig_ref = path[path.CMapId == qrycontigid]
    else:
        logger.error(
            "Something went wrong in reading molecules, recheck BNGO folder")
        raise KeyError(
            "Something went wrong in reading molecules, recheck BNGO folder")

    if mult == -1:
        # SiteID always contains one more label which is the end of contig
        last_label = _cmap_contig_ref["SiteID"].max()-1
        return last_label - last_label_query
    logger.error(
        "Something went wrong in reading molecules, recheck BNGO folder")
    raise KeyError(
        "Something went wrong in reading molecules, recheck BNGO folder")


def _gnoucl(subrow: pd.Series, molecules: pd.DataFrame, telarm: Literal["left", "right"],
            master_orientation=Literal["+", "-"],
            chrom_query: str = CHROM_QUERY) -> pd.Series:
    """
    A wrapper around 'get_number_of_unpaired_contig_labels'

    """
    aln_str = subrow["RefContigMolPair"]
    qry_id = subrow["QryContigID"]
    orient = subrow["Orientation"]
    return get_number_of_unpaired_contig_labels(path=molecules,
                                                alignment=aln_str,
                                                qrycontigid=qry_id,
                                                orientation=orient,
                                                telarm=telarm,
                                                master_orientation=master_orientation,
                                                chrom_query=chrom_query)
# %%


def alignments_to_reference(xmap_df: pd.DataFrame, label: int) -> pd.DataFrame:
    """
    An XMAP dataframe contains column Alignment which is a list
    of tuples of the type:
    >>>> (ref_label_A,query_label_A)(ref_label_B,query_label_B)
    Which is pair of aligned labels.

    Given a label, we want to find only those parts of a XMAP dataframe
    that contain mapping onto said label as a reference
    """
    return xmap_df.loc[xmap_df.Alignment.str.contains(f"({label},", regex=False)]


def extract_reference_query_pair(xmap_df: pd.DataFrame, aligned_label: int) -> pd.Series:
    """
    An XMAP dataframe contains column Alignment which is a list
    of tuples of the type:
    >>> (ref_label_A,query_label_A)(ref_label_B,query_label_B)
    Which is pair of aligned labels.

    Given a label onto the reference, we want to find its pair on
    the query.

    For example, if we have an Alignment of:
    >>>> (1,15)(2,16)(3,17)(4,15)
    and aligned_label of 3, the pair is 17

    """

    return xmap_df.Alignment.str.extract(r"\((%d),(\d+)\)" % aligned_label)[1].astype(int)


def read_contig_xmap(contig_path: str, aligned_label: Optional[int] = None) -> pd.DataFrame:
    """
    Given an XMAP file found at contig_path,
    this function reads the XMAP file and given a label position on the
    reference (for example 2),
    this function then returns only parts of XMAP files
    that are aligned to said label.
    """
    # read contig in the XMAP format
    contig_aligned = read_map_file(contig_path)
    if aligned_label is not None:
        # find only alignments to our first_aligned_query
        contig_aligned = alignments_to_reference(contig_aligned, aligned_label)
    return contig_aligned


def read_molecules_cmap(molecules_path: str, molecule_ids: Optional[Iterable[int]] = None) -> pd.DataFrame:
    """
    Given an CMAP file found at molecules path,
    this function reads the CMAP file and then
    if an iterable is given, returns parts of those CMAP
    file that contain given molecules.

    Note: It just returns "CMapId","SiteID", and "Position" columns
    """

    molecules = read_map_file(molecules_path)
    if molecule_ids is not None:
        # E1101 is a false positive
        molecules = molecules.loc[molecules.CMapId.isin(
            molecule_ids)]  # pylint:disable=E1101
    return molecules[["CMapId", "SiteID", "Position"]]


def remap_query_position(contig_aligned: pd.DataFrame, molecules: pd.DataFrame, aligned_label: int) -> pd.DataFrame:
    """
    Given an XMAP in contig_aligned, a CMAP in molecules, and aligned label
    this function first extracts a  query paired to said reference label and then
    finds nucleotide position of said query


    For example, suppose that you have an alignment of:
    >>> (1,15)(2,16)(3,17)(4,15)
    In that case, QryStartPos will be the position of label 15, as the label
    15 is paired to the label 1, and label 1 is the first pairing, for example
    >>> 15 345.2

    But if you give aligned_label=3, then it will find that the query pair label
    is 17, and it will update QryStartPos to be the position of the 17th label,
    >>> 25 343.2

    The new Positions are stored in the column called "Position" so
    that the function works with both contig-to-reference orientations:
        for + it's QryEndPos
        for - it's QryStartPos
    """
    # finds the query pair of given reference label
    contig_aligned["LastLabelPair"] = extract_reference_query_pair(
        contig_aligned, aligned_label)

    # merges information on from molecules
    # drops extraneous info - CMapId, SiteID and LastLabelPair
    contig_aligned = pd.merge(contig_aligned, molecules,
                              left_on=["QryContigID", "LastLabelPair"],
                              right_on=["CMapId", "SiteID"]
                              ).drop(["CMapId", "SiteID", "LastLabelPair"], axis=1)
    # I CAN DROP QRY START POS AND RETURN POSITION
    # FOR INTEGRATION BETWEEN THE NEGATIVE AND POSITIVE
    return contig_aligned


# %%
def calculate_telomere_length_formula(row: pd.Series, telarm: Literal["left", "right"],
                                      master_orientation: Literal["+", "-"]) -> pd.Series:
    """
    This function calculates the length of telomere.

    TODO - better documentation
    """

    orientation = row["Orientation"]
    qry_len = row["QryLen"]
    pos = row["Position"]

    map_orient = {"+": 1, "-": -1}[orientation]
    map_telarm = {"left": 1, "right": -1}[telarm]
    map_morient = {"+": 1, "-": -1}[master_orientation]

    mult = map_orient * map_telarm * map_morient

    if mult == 1:
        return pos
    if mult == -1:
        return qry_len - pos
    logger.error("Something went wrong in orientation")
    raise KeyError("Something went wrong in orientation")


def calculate_telomere(row: pd.Series, path: str,
                       contig_xmap_format: str = CONTIG_XMAP,
                       contig_query_format: str = CONTIG_QUERY,
                       chrom_query: str = CHROM_QUERY,):

    telarm: Literal["left", "right"] = row["TelomereArm"]
    # master_orientation:Literal["+","-"] = row["Orientation"]
    # master_chromosome = row["RefContigID"]
    master_contig: int = row["QryContigID"]
    contig_alignment: str = row["Alignment"]

    # pathings of contig XMAP and contig as _q.cmap
    contig_path: str = joinpaths(
        path, contig_xmap_format.format(x=master_contig))
    molecules_path: str = joinpaths(
        path, contig_query_format.format(x=master_contig))

    # get the label of bount reference
    aligned_label: int = get_first_query_bound_reference(
        contig_alignment, telarm)

    # reads only those alignments that contain query label
    # on the reference
    contig_aligned: pd.DataFrame = read_contig_xmap(contig_path, aligned_label)

    # reads molecules found in alignment to the last reference
    molecules: pd.DataFrame = read_molecules_cmap(
        molecules_path, contig_aligned.QryContigID.unique())

    # Given a label onto the reference, we want to find its pair on the query.
    # For example, if we have an Alignment of:
    # >>>> (1,15)(2,16)(3,17)(4,15)
    # and aligned_label of 3, the pair is 17
    contig_aligned["RefContigMolPair"] = contig_aligned.Alignment.str.extract(
        r"\((%d),(\d+)\)" % aligned_label)[1].astype(str)

    # we then format it as a pair (3,17)
    contig_aligned["RefContigMolPair"] = f"({aligned_label}," + \
        contig_aligned["RefContigMolPair"] + ")"

    # we extract unpaired labels
    contig_aligned["UnpairedReferenceLabels"] = row["OffsetLabel"]
    contig_aligned["UnpairedContigLabels"] = get_number_of_unpaired_contig_labels(path=path,
                                                                                  alignment=row["Alignment"],
                                                                                  qrycontigid=row["QryContigID"],
                                                                                  orientation=row["Orientation"],
                                                                                  chrom_query=chrom_query, telarm=telarm,
                                                                                  master_orientation="+",  # necessary for the formula
                                                                                  )

    contig_aligned["UnpairedMoleculeLabels"] = contig_aligned.apply(_gnoucl, axis=1,   # pylint:disable=E1101,E1137
                                                                    molecules=molecules, telarm=telarm,
                                                                    master_orientation=row["Orientation"])

    # calculate telomeres
    # first remap position
    contig_aligned = remap_query_position(
        contig_aligned, molecules, aligned_label)

    # then perform formula
    contig_aligned["TelomereLen"] = contig_aligned.apply(
        calculate_telomere_length_formula, axis=1, telarm=telarm,
        master_orientation=row["Orientation"])

    # refactor contig for easier data insertion
    contig_aligned = contig_aligned[["QryContigID", "RefContigID", "Orientation", "Confidence",
                                     "QryLen", "TelomereLen",
                                     "UnpairedMoleculeLabels", "UnpairedContigLabels", "UnpairedReferenceLabels"]
                                    ]

    contig_aligned.columns = ["MoleculeID", "QryContigID", "MoleculeOrientation", "MoleculeConfidence",
                              "MoleculeLen", "TelomereLen",
                              "UnpairedMoleculeLabels", "UnpairedContigLabels", "UnpairedReferenceLabels"]
    # inserting information
    contig_aligned["RefContigID"] = row["RefContigID"]
    contig_aligned["ContigOrientation"] = row["Orientation"]
    contig_aligned["BoundReferencePosition"] = row["BoundReferencePosition"]

    if telarm == "left":
        contig_aligned["AlignedLabelPosition"] = row["RefStartPos"]
        contig_aligned["EndDistance"] = row["RefStartPos"]

    if telarm == "right":
        contig_aligned["AlignedLabelPosition"] = row["RefEndPos"]
        contig_aligned["EndDistance"] = row["RefLen"] - row["RefEndPos"]

    contig_aligned["Telomere"] = telarm
    contig_aligned["ContigOrientation"] = contig_aligned["ContigOrientation"].replace({
                                                                                      "+": 1, "-": -1})
    contig_aligned["MoleculeOrientation"] = contig_aligned["MoleculeOrientation"].replace({
                                                                                          "+": 1, "-": -1})
    # reorder columns
    columns = ["RefContigID", "MoleculeID",  "QryContigID",
               "MoleculeOrientation", "ContigOrientation", "MoleculeConfidence",
               "AlignedLabelPosition", "BoundReferencePosition",
               "UnpairedReferenceLabels", "UnpairedContigLabels", "UnpairedMoleculeLabels",
               "EndDistance",
               "MoleculeLen", "Telomere", "TelomereLen",
               ]
    return contig_aligned[columns]


def time_function(*args, **kwargs) -> Callable:
    "timing telomere calculation, adjusting decorating for multiprocessing"
    return func_timer(calculate_telomere)(*args, **kwargs)

# %%


def _calculate_telomere_lengths(path: str, how: Literal["left", "right"],
                                chrom_xmap: str = CHROM_XMAP, chrom_reference: str = CHROM_REFERENCE,
                                contig_xmap_format: str = CONTIG_XMAP, contig_query_format: str = CONTIG_QUERY,
                                chrom_query: str = CHROM_QUERY, threads: int = 1,
                                ) -> pd.DataFrame:
    """
    Calculates telomere lengths for a given path to BNGO de novo Assembly and
    returns list of DataFrames

    """
    # TODO - better documentation
    chrom_xmapdf = fish_paired_label(path=path, how=how,
                                     chrom_xmap=chrom_xmap, chrom_reference=chrom_reference)

    # chrom_xmapdf = chrom_xmapdf[chrom_xmapdf.RefContigID.isin([4,14,20,21])]
    # chrom_xmapdf = chrom_xmapdf.head(1)

    frozen_calculation = partial(time_function, path=path, contig_xmap_format=contig_xmap_format,
                                 contig_query_format=contig_query_format, chrom_query=chrom_query)
    iterator = (x[1] for x in chrom_xmapdf.iterrows())

    if threads > 1:
        with Pool(threads) as pool:
            results = pool.map(frozen_calculation, iterator)
    else:
        results = [frozen_calculation(row) for row in iterator]

    return results


def reduce_dataset(data: pd.DataFrame, ref_tol: int,
                   con_tol: int, mol_tol: int, dis_tol: int) -> pd.DataFrame:
    """
    Reduces our dataset based on the numbers provided.
    For every alignment we have the following three statistics:
        UnpairedReferenceLabels
        UnpairedContigLabels
        UnpairedMoleculeLabels

    Also, additional parameter ('dis_tol') checks 'EndDistance' which
    is the distance between first aligned label and the end of the chromosome.
    In case where our label is too far from the chromosome ends, we can simply
    exclude those molecules from calculation.

    Our alignment is MOLECULE to CONTIG to REFERENCE and this alignment
    doesn't have to be perfect. For example, we can align the penultimate
    label of the CONTIG to last label of REFERENCE.

    These three statistics tell us how many labels are after the last aligned
    pair. Depending on your level of tolerance, this can vary.

    For example, passing 'ref_tol'=1 means that you tolerate one extra label
    after the last paired label on the reference.
    See recommended values in .const.
    """
    input_df = data.copy()

    original_len = len(data)
    len_input_df = original_len

    cols = ["EndDistance", "UnpairedReferenceLabels",
            "UnpairedContigLabels", "UnpairedMoleculeLabels"]
    for name, tol in zip(cols, [dis_tol, ref_tol, con_tol, mol_tol]):
        input_df = input_df.loc[input_df[name] <= tol]
        len_new = len(input_df)

        logger.info("Excluded %s > %d from input_df - went from %d to %d molecules (%.2f percent reduction)",
                    name, tol, len_input_df, len_new, 100 - len_new/len_input_df*100)
        len_input_df = len_new

        if len_new == 0:
            break

    logger.info("Went from %d to %d molecules (%.2f percent reduction)",
                original_len, len_new, 100 - len_new/original_len*100
                )
    return input_df
# %% MOLECULAR STATISTICS


def correlation_test(telomeres: pd.DataFrame, name: str = "") -> None:
    "tests pearson correlation between length of molecule and telomere"
    corr, pval = stats.pearsonr(
        telomeres["MoleculeLen"], telomeres["TelomereLen"])
    pval = np.log10(pval)
    logger.info(
        "%s - LABEL_STATS - MoleculeLen/TelomereLen correlation %.4f with log(pvalue) %.8f ", name, corr, pval)


def get_average_mol_len(path: str, jason: str = INFREP) -> Optional[float]:
    "reads average molecule length from informatics report"
    try:
        with open(joinpaths(path, jason), "r",  encoding="utf-8") as _f:
            # Load the JSON data into a Python dictionary
            data = json.load(_f)
    except FileNotFoundError:
        logger.error("No informatics report found")
        return None

    inp_mol_stats = data.get("Input molecule stats (filtered)", None)

    if not isinstance(inp_mol_stats, dict):
        logger.error("Key %s not found", "Input molecule stats (filtered)")
        return None

    av_len = inp_mol_stats.get("Average length (kbp)", None)

    if av_len is None:
        logger.error("Key %s not found", "Average length (kbp)")
        return None

    return float(av_len) * 1000  # control for kilobases


def t_test_mols(telomeres: pd.DataFrame, path: str, name: str = "",
                jason: str = INFREP) -> None:
    """
    tests if telomeric molecules are different in length from the
    average molecule
    """
    avlen = get_average_mol_len(path, jason)
    if avlen is None:
        return None

    t_stat, pval = stats.ttest_1samp(telomeres["MoleculeLen"], avlen)

    # one-sided cohen_d is just dividing t_stat by square root of N
    cohen_d = t_stat / len(telomeres)**0.5

    pval = np.log10(pval)
    logger.info("%s - LABEL_STATS - Difference from expected molecule length -  Cohen D is %.4f with log(pvalue) %.8f ", name, cohen_d, pval)


def molecule_statistics(telomeres: pd.DataFrame, path: str, name: str = "",
                        jason: str = INFREP) -> None:
    """
    Calculates relevant statistics:
        1. Correlation of TelomereLen with MoleculeLen
        2. P value of MoleculeLen with expected value of molecule len being
        in informatics report along with CohenD
    """
    correlation_test(telomeres, name)
    t_test_mols(telomeres, path, name, jason)


def calculate_telomere_lengths(path: str,
                               how: Union[List[Literal["left"]], List[Literal["right"]], List[Literal["left", "right"]]],
                               chrom_xmap: str = CHROM_XMAP, chrom_reference: str = CHROM_REFERENCE,
                               contig_xmap_format: str = CONTIG_XMAP, contig_query_format: str = CONTIG_QUERY,
                               chrom_query: str = CHROM_QUERY, threads: int = 1,
                               ref_tol: int = REF_TOL, con_tol: int = CON_TOL, mol_tol: int = MOL_TOL,
                               dis_tol: int = DIS_TOL
                               ) -> pd.DataFrame:

    logger.info("Calculating telomere length for file found at %s", path)

    final_result = list()
    name = os.path.basename(path)

    for option in how:

        logger.info("Calculating telomere lengths for %s arm ", option.upper())

        data = _calculate_telomere_lengths(path=path, how=option,
                                           chrom_xmap=chrom_xmap, chrom_reference=chrom_reference,
                                           contig_xmap_format=contig_xmap_format, contig_query_format=contig_query_format,
                                           chrom_query=chrom_query, threads=threads,
                                           )

        data = pd.concat(data, axis=0, ignore_index=True)

        # display statistics about unpaired labels
        _cols = ["RefContigID", "QryContigID", "UnpairedReferenceLabels", "ContigOrientation",
                 "UnpairedContigLabels", "EndDistance"]
        stats_df = data[_cols].drop_duplicates()

        for _, row in stats_df.iterrows():
            display = ' , '.join(
                [f"{index} : {value:.0f}" for index, value in row.items()])
            logger.info("%s - LABEL_STATS - %s", name, display)
        filtered_dataset = reduce_dataset(
            data, ref_tol=ref_tol, con_tol=con_tol, mol_tol=mol_tol, dis_tol=dis_tol)
        final_result.append(filtered_dataset)

    if len(final_result) > 1:
        final_result = pd.concat(final_result, ignore_index=True)
    else:
        final_result = final_result[0]
    # perform final testing
    molecule_statistics(final_result, path, name, INFREP)
    return final_result
