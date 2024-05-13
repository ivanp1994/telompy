# -*- coding: utf-8 -*-
"""
Created on Wed May  8 10:11:51 2024

@author: ivanp
"""
from typing import Iterable, Callable, Union, Literal, Optional
from functools import partial
from multiprocessing import Pool

import pandas as pd

from .utils import read_map_file, joinpaths, func_timer
from .const import CONTIG_PATH, QUERYCMAP_PATH, MASTER_XMAP, MASTER_REFERENCE, MASTER_QUERY
from .const import LOGGER as logger


pd.options.mode.chained_assignment = None  # default='warn'

# %% FUNCTIONS FOR FISHING LAST LABELS ON MASTER ASSEMBLY


def fish_last_aligned_label(xmap_df: pd.DataFrame) -> pd.DataFrame:
    """
    From a given xmap_df we find the last aligned reference label
    for every alignment - and then we filter out contigs ('QryContigID')
    that do not align to said last aligned reference.

    If there are N chromosomes (referent contigs - 'RefContigID'), there
    should be M (M>=N) rows in the resulting output
    """
    last_aligned = xmap_df.copy()

    # last aliged labels
    # pattern = r"\((\d+),(\d+)\)(?!.*\(\d+,\d+\))"
    # master_xmap["LastLabel"] = master_xmap["Alignment"].str.extract(pattern)[0].astype(int)
    # THE ABOVE IS TOO SLOW

    last_aligned["SiteID"] = last_aligned["Alignment"].str.split(")(", regex=False).str[-1]  # finds last element
    last_aligned["SiteID"] = last_aligned["SiteID"].str.split(",").str[0].astype(int)  # finds first pair of last element

    # last_aligned.drop("Alignment", axis=1, inplace=True)
    # above is commented 'cause alignment is horrible to parse in pickle

    # finds the last last aligned pair for every contig
    last_aligned = last_aligned.groupby("RefContigID").apply(lambda one_chrom: one_chrom[one_chrom["SiteID"] == one_chrom["SiteID"].max()])
    last_aligned = last_aligned.reset_index(drop=True)
    return last_aligned


def fish_offset(last_aligned: pd.DataFrame, reference_cmap: pd.DataFrame) -> pd.DataFrame:
    """
    Given an XMAP file in last_aligned that only has query contigs aligned
    to the last aligned reference label and a CMAP file in reference_cmap this
    function finds the difference between the last label on reference and
    the last aligned label on the reference.

    For example, suppose a QryContigID of 15 is aligned to RefContigID of 3
    and the last alignment is (26700,2) - meaning that the label 26700
    on the reference contig is aligned to label 2 on the query contig.

    The position of 26700th label is 18,365,800 as found in reference_cmap.
    If 26700th label is truly the last label on the reference contig, the
    'Offset' value will be 0.

    But if we have 267001st label whose position is 19,365,800, meaning that
    the last aligned label on reference is not the last label on reference, then
    said offset will be 19,365,800 - 18,365,800 = 1,000,000

    EXTRA: added info about the number of Labels on Reference not mapped to
    contig

    """
    # last_aligned = fish_last_aligned_label(master_xmap)
    # get positions of last aligned labels
    aligned_positions = last_aligned[["RefContigID", "SiteID", "QryContigID"]].drop_duplicates()

    aligned_positions = pd.merge(left=aligned_positions,
                                 right=reference_cmap[["CMapId", "SiteID", "Position"]],
                                 left_on=["RefContigID", "SiteID"],
                                 right_on=["CMapId", "SiteID"],
                                 )[["RefContigID", "QryContigID", "SiteID", "Position"]]

    # get positions of last labels
    # last label on reference
    last_labels = reference_cmap[reference_cmap.SiteID == reference_cmap.NumSites][["CMapId", "SiteID", "Position"]]
    last_labels.columns = ["RefContigID", "LastID", "AlignedLabelPosition"]

    # merge the two
    aligned_positions = pd.merge(left=aligned_positions,
                                 right=last_labels,
                                 left_on="RefContigID",
                                 right_on="RefContigID",
                                 how="left")

    #aligned_positions["Offset"] = aligned_positions["LastPosition"] - aligned_positions["Position"]
    aligned_positions["OffsetLabel"] = aligned_positions["LastID"] - aligned_positions["SiteID"]
    return aligned_positions


def fish_last_label(path: str, main_xmap: str = MASTER_XMAP, main_cmapr: str = MASTER_REFERENCE) -> pd.DataFrame:
    """
    Given a path to the BioNano denovo assembly, this function finds XMAP
    for assembly of contigs-to-chromosome, then finds alignments that
    align to the last (the largest label on the reference) label on the reference.

    In a situation where the last aligned label on the reference is not the
    last label on reference, an additional column called 'Offset' is found.
    This contains the difference between last label and last aligned label
    in base pairs.
    """

    master_xmap = read_map_file(joinpaths(path, main_xmap))
    master_cmap = read_map_file(joinpaths(path, main_cmapr))

    # find last aligned reference label on the master xmap
    # and filter out only contigs that align to it
    last_aligned = fish_last_aligned_label(master_xmap)

    # find the distance between the last aligned reference label
    # and the last reference label
    # if a contig maps onto the last label then the OffSet is 0
    aligned_offsets = fish_offset(last_aligned, master_cmap)

    # merge them
    last_aligned = pd.merge(left=last_aligned,
                            #right=aligned_offsets[["QryContigID", "RefContigID", "Offset", "Offset_Label","AlignedLabelPosition"]],
                            right=aligned_offsets[["QryContigID", "RefContigID", "OffsetLabel","AlignedLabelPosition"]],
                            left_on=["QryContigID", "RefContigID"],
                            right_on=["QryContigID", "RefContigID"])
    return last_aligned

# %% FUNCTIONS FOR CALCULATING TELOMERE LENGTHS


def get_first_query_last_reference(align_str: str) -> int:
    """
    An XMAP dataframe contains column Alignment which is a list
    of tuples of the type:
    >>>> (ref_label_A,query_label_A)(ref_label_B,query_label_B)
    Which is pair of aligned labels.

    We want to find the first aligned label on Query that
    aligns on the last aligned label on Reference for cases when
    the first label (as in label number 1) does not actually map
    onto the last label on the reference.

    For that purpose, we convert the string series to a tuple of integers,
    and then we use Python's convention of sorting tuples:
    when finding maximum value of a list of tuples, only
    first the value is considered

    So we return the Query label of last Reference label

    What if the last label aligns to multiple labels???


    The XMAP documentation says:
        : When two sites in the reference align with the same site in the query
        , it is an indication that the two sites in the reference failed to
        resolve. Alignment provides a view of aligned pairs which would normally
        be ignored by HitEnum (CIGAR string).

    So even if the sites in reference failed to resolve, I guess it wont matter
    much for us.
    """
    # ignoring R1728
    # https://github.com/pylint-dev/pylint/pull/3309#discussion_r576683109
    # generator is slower for smaller lists - and we have roughly maximum of ~20k
    # microbenchmarking I've gotten generator to be ~3 ms slower
    align = [tuple([int(x) for x in x.split(",")]) for x in align_str[1:-1].split(")(")]  # pylint: disable=R1728
    return max(align)[1]


def get_number_of_unpaired_contig_labels(path: Union[str, pd.DataFrame], alignment: str, qrycontigid: int,
                                         orientation: Literal["-", "+"],
                                         contig_query: str = MASTER_QUERY) -> int:
    """
    This function returns number of unpaired labels on the contig (query).
    The contig is a query in a contig-to-reference assembly.

    For example, suppose you have a pairing:
        (.,.)(.,.)...(27777,35)
    The label 27777 on the reference is paired with label 35 on the query.

    How many labels on the query are unpaired with the reference?
    If the orientation is negative, then it means that the labels are in a
    descending order:
        ...(27775,37),(27776,36),(27777,35)
    Meaning that there are 34 labels that are unpaired.

    However, if the orientation is positive, then it means that the labels
    are in an ascending order:
        ... (27775,33),(27776,34),(27777,35)

    We cannot know how many (if any) labels are unpaired.
    In that case, we must open the query CMAP, and find the last label.
    For example, if we found that the last label is 50,
    then it means we have 15 unpaired labels.



    NOTE: It returns number of labels AFTER last pairing, not
    TOTAL number of unpaired labels. For example:
        (1,4)(2,2) of - Orientation
    Contains only one unpaired label after the last pair, and that is label 1.
    Label 3 is also unpaired, but not taken into consideration for the
    purpose of this project
    """
    if orientation not in ["-", "+"]:
        raise ValueError("Invalid orientation - must be either '-' or '+'")

    last_label_qry = get_first_query_last_reference(alignment)

    # if orientation is negative
    # then queries are descending
    # so it always ends with 1
    if orientation == "-":
        return last_label_qry - 1

    # if it's positive, we need to info about molecules stored in query
    # if the path is a string
    if isinstance(path, str):
        _cmap_contig_ref = read_molecules_cmap(joinpaths(path, contig_query), [qrycontigid])

    # if the path is an already loaded CMAP
    elif isinstance(path, pd.DataFrame):
        _cmap_contig_ref = path[path.CMapId == qrycontigid]

    # raise ValueError otherwise
    # TODO - maybe refactor this into logging?
    else:
        raise ValueError("Path must be either a string or a DataFrame")

    last_label = _cmap_contig_ref["SiteID"].max()-1
    return last_label - last_label_qry


def _gnoucl(subrow: pd.Series, molecules: pd.DataFrame) -> pd.Series:
    """
    A wrapper around 'get_number_of_unpaired_contig_labels'

    """
    aln_str = subrow["RefContigMolPair"]
    qry_id = subrow["QryContigID"]
    orient = subrow["Orientation"]
    return get_number_of_unpaired_contig_labels(path=molecules,
                                                alignment=aln_str,
                                                qrycontigid=qry_id,
                                                orientation=orient)


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
        molecules = molecules.loc[molecules.CMapId.isin(molecule_ids)]  # pylint:disable=E1101
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
    contig_aligned["LastLabelPair"] = extract_reference_query_pair(contig_aligned, aligned_label)

    # merges information on from molecules
    # drops extraneous info - CMapId, SiteID and LastLabelPair
    contig_aligned = pd.merge(contig_aligned, molecules,
                              left_on=["QryContigID", "LastLabelPair"],
                              right_on=["CMapId", "SiteID"]
                              ).drop(["CMapId", "SiteID", "LastLabelPair"], axis=1)
    # I CAN DROP QRY START POS AND RETURN POSITION
    # FOR INTEGRATION BETWEEN THE NEGATIVE AND POSITIVE
    return contig_aligned


def calculate_telomere(row: pd.Series, path: str,
                       contig_format: str = CONTIG_PATH,
                       querycmap_format: str = QUERYCMAP_PATH,
                       contig_query: str = MASTER_QUERY,
                       gap_size:Optional[int]=None) -> pd.DataFrame:
    """
    Given a row, which is a pandas Series, and path which is the location
    of a BNGO assembly, calculates telomere length.

    The output is a pandas DataFrame with columns:
        RefContigID:
            The 'chromosome' (highest level assembly)
        MoleculeID:
            The ID of a molecule that maps to a contig
        QryContigID:
            The ID of a contig that consists of molecules and maps
            to the RefContigID
        MoleculeConfidence:
            The Confidence of assembly of molecule-to-contig
        TelomereLen:
            Length of telomere as calculated by the formula
        TelomereLen_corr:
            Length of telomere reduced by the distance of the last ALIGNED label
            on the reference from the last label on the reference
            When the last aligned label on the reference is equal to the last
            aligned label on the refrence, this is equal to TelomereLen


    """
    # TODO -- add formula
    # TODO - better documentation
    # ROW UNPACKING
    master_orientation = row["Orientation"]
    # master_chromosome = row["RefContigID"]
    master_contig = row["QryContigID"]
    contig_alignment = row["Alignment"]

    logger.info("Processing telomeres mapped to RefContigID %d from QryContigID %d",
                row["RefContigID"], row["QryContigID"])

    # READING CONTIG PATHS
    contig_path = joinpaths(path, contig_format.format(x=master_contig))
    molecules_path = joinpaths(path, querycmap_format.format(x=master_contig))

    # WORKFLOW
    # gets first label on query that is mapped to the last reference
    aligned_label = get_first_query_last_reference(contig_alignment)
    # reads only those alignments that contain query label
    # on the reference
    contig_aligned = read_contig_xmap(contig_path, aligned_label)

    # reads molecules found in alignment to the last reference
    molecules = read_molecules_cmap(molecules_path, contig_aligned.QryContigID.unique())

    # returns number of labels on molecules that are unpaired
    paired_query = extract_reference_query_pair(contig_aligned, aligned_label)
    contig_aligned["RefContigMolPair"] = paired_query.astype(str).apply(lambda x: f"({aligned_label},{x})")  # pylint:disable=E1137
    contig_aligned["UnpairedMoleculeLabels"] = contig_aligned.apply(_gnoucl, axis=1,   # pylint:disable=E1101,E1137
                                                                    molecules=molecules)
    # pulling out position of pair of the reference label
    contig_aligned = remap_query_position(contig_aligned, molecules, aligned_label)

    # extracting length of telomeres according to the formula
    contig_aligned["TelomereLen"] = contig_aligned["Position"]
    contig_aligned.loc[contig_aligned["Orientation"] == master_orientation,
                       "TelomereLen"] = contig_aligned["QryLen"] - contig_aligned["Position"]  # pylint:disable=C0301
    # correct for the offset
    #contig_aligned["TelomereLen_corr"] = contig_aligned["TelomereLen"] - row["Offset"]

    # reformat and return
    # NOTE: deleted TelomereLenCorr
    contig_aligned = contig_aligned[["QryContigID", "RefContigID", "Orientation",
                                     "Confidence", "TelomereLen",  "UnpairedMoleculeLabels","QryLen"]]
    contig_aligned.columns = ["MoleculeID", "QryContigID", "MoleculeOrientation",
                              "MoleculeConfidence", "TelomereLen",  "UnpairedMoleculeLabels","MoleculeLen"]

    # inserting information
    contig_aligned.insert(0, "RefContigID", row["RefContigID"])
    contig_aligned.insert(4, "ContigOrientation", master_orientation)
    contig_aligned.insert(5,"AlignedLabelPosition",row["AlignedLabelPosition"])
    contig_aligned.insert(6,"LastReferencePosition",row["RefEndPos"])

    contig_aligned["UnpairedReferenceLabels"] = row["OffsetLabel"]
    contig_aligned["UnpairedContigLabels"] = get_number_of_unpaired_contig_labels(path=path,
                                                                                  alignment=row["Alignment"],
                                                                                  qrycontigid=row["QryContigID"],
                                                                                  orientation=row["Orientation"],
                                                                                  contig_query=contig_query
                                                                                  )

    
    # telomere correction via gap size
    
    if gap_size is None:
        gap_size = row["RefLen"] - row["AlignedLabelPosition"]
    #TODO - do I insert last aligned label position
    offset = row["AlignedLabelPosition"] - row["RefLen"] + gap_size
    contig_aligned["TelomereLen_corr"] = contig_aligned["TelomereLen"] + offset

    # a vain attempt to decrease memory usage by casting length to an integer
    contig_aligned["TelomereLen"] = contig_aligned["TelomereLen"].astype(int)
    contig_aligned["TelomereLen_corr"] = contig_aligned["TelomereLen_corr"].astype(int)

    return contig_aligned


def time_function(*args, **kwargs) -> Callable:
    "timing telomere calculation, adjusting decorating for multiprocessing"
    return func_timer(calculate_telomere)(*args, **kwargs)


def calculate_telomere_lengths(path: str, main_xmap: str = MASTER_XMAP, main_cmapr: str = MASTER_REFERENCE,
                               contig_format: str = CONTIG_PATH, querycmap_format: str = QUERYCMAP_PATH,
                               threads: int = 1, gap_size:Optional[int]=None) -> Iterable[pd.DataFrame]:
    """
    Calculates telomere lengths for a given path to BNGO de novo Assembly and
    returns list of DataFrames

    """
    # TODO - better documentation
    main_xmapdf = fish_last_label(path, main_xmap, main_cmapr)

    # main_xmapdf = main_xmapdf[main_xmapdf.RefContigID.isin([4,14,20,21])]
    # main_xmapdf = main_xmapdf.head(1)

    frozen_calculation = partial(time_function, path=path, contig_format=contig_format,
                                 querycmap_format=querycmap_format,gap_size=gap_size)
    iterator = (x[1] for x in main_xmapdf.iterrows())

    if threads > 1:
        with Pool(threads) as pool:
            results = pool.map(frozen_calculation, iterator)
    else:
        results = [frozen_calculation(row) for row in iterator]
    return results
    
    #todo - pass concat into outer function
    # try:
    #     concatenated = pd.concat(results, ignore_index=True)
    # except MemoryError:
    #     concatenated = results
    # return concatenated
    
### THE PROBLEM WITH THE ABOVE IS THAT IT PRODUCES
