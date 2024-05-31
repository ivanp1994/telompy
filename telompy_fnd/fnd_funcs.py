# -*- coding: utf-8 -*-
"""

@author: ivanp
"""
import os

import logging
from typing import (List, Set, Generator,
                    Tuple, Dict,
                    Optional, Literal,Union)

import numpy as np
import pandas as pd

from .utils import (func_timer,  get_partial_alignment,
                    count_comment_lines, read_map_file)
from .const import REF_TOL,MOL_TOL,DIS_TOL

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)
logger = logging.getLogger("telompy_fandom")

TelarmType = Union[
    List[Literal["left"]],
    List[Literal["right"]],
    List[Literal["left", "right"]]
]
#%% FUNCTIONS REUSED FROM BNGO BRANCH

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
    
    
    return aligned_positions # it returns Position


def fish_paired_label(master_xmap: pd.DataFrame, master_cmap: pd.DataFrame, how: Literal["left", "right"],
                      ) -> pd.DataFrame:
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

#%% FUNCTIONS FOR FANDOM BRANCH

@func_timer
def read_fandom_xmap(xmap_path: str, ref_path: str,
                     how: TelarmType = ["left", "right"],
                     chunksize: int = 500000, ) -> pd.DataFrame:

    master_cmap = read_map_file(ref_path)
    skiprows = count_comment_lines(xmap_path)
    elements = list()
    
    with pd.read_csv(xmap_path, chunksize=chunksize, sep="\t",
                     skiprows=skiprows, usecols=[1, 2,5,6, 7,8, 11,13,],
                     names=["QryContigID", "RefContigID", "RefStartPos","RefEndPos",
                            "Orientation", "Confidence","RefLen","Alignment"]
                     ) as reader:

        for chunk in reader:
            if "left" in how:
                elements.append(fish_paired_label(chunk, master_cmap, "left"))
            if "right" in how:
                elements.append(fish_paired_label(chunk, master_cmap, "right"))
    data = pd.concat(elements)

    # this corrects for chunksize reading on XMAP
    # what can happen is that in one chunk, no chrom 1 has 3 offset label
    # and in another chunk the chrom 1 has 1 offset label
    # we want to keep chrom 1 from the second chunk

    data = data.groupby(["RefContigID", "TelomereArm"]).apply(
        lambda one_chrom: one_chrom[one_chrom["OffsetLabel"] == one_chrom["OffsetLabel"].min()])

    data["Orientation"] = data["Orientation"].map({"+": 1, "-": -1})

    data["TelomereArm"] = data["TelomereArm"].map({"left": 1, "right": -1})


    # TODO - maybe do it smarter, without str
    # but would apply be slower?
    # does it matter?
    # find the pair of the SiteID column and put it next to it
    # Use regex to extract all pairs into lists of tuples
    data['AlignmentPairs'] = data['Alignment'].str.findall(r'\((\d+),(\d+)\)')
    # Convert the lists of tuples from strings to integers
    data['AlignmentPairs'] = data['AlignmentPairs'].apply(lambda x: [(int(a), int(b)) for a, b in x])
    # Create the "QueryLabel" column by finding the corresponding second element for each SiteID
    data['QueryLabel'] = data.apply(lambda row: next((b for a, b in row['AlignmentPairs'] if a == row['SiteID']), None), axis=1)
    # Drop the intermediate 'AlignmentPairs' column
    data = data.drop(columns=['AlignmentPairs'],)
    
    # now calculate the distance from the chromosome' ends to first/last aligned label
    data['EndDistance'] = np.where(
    data['TelomereArm'] == 1,
    data['RefStartPos'],
    data['RefLen'] - data['RefEndPos']
    )
    
    # also get AlignedLabelPosition - position of first/last aligned label on the chromosome
    data['AlignedLabelPosition'] = np.where(
    data['TelomereArm'] == 1,
    data['RefStartPos'],
    data['RefEndPos']
    )
    # select only specific columns
    cols = ["QryContigID", "RefContigID", "Orientation","Confidence","TelomereArm","OffsetLabel", "QueryLabel","EndDistance","AlignedLabelPosition","BoundReferencePosition"]
    data = data[cols]
    # rename those columns
    cols = ["QryContigID", "RefContigID", "MoleculeOrientation","MoleculeConfidence","TelomereArm","UnpairedReferenceLabels", "QueryLabel","EndDistance","AlignedLabelPosition","BoundReferencePosition"]
    data.columns = cols
   

    return data[cols].reset_index(drop=True)


# %% bnx function

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

#%% telomere function
@func_timer
def calculate_telomeres(telomeric_mols:List[np.array],telomeric_align_df:pd.DataFrame)->pd.DataFrame:
    
    # get columns relevant for calculation
    qry_id_col = telomeric_align_df.columns.get_loc("QryContigID")
    orientation_col = telomeric_align_df.columns.get_loc("MoleculeOrientation")
    telomere_col = telomeric_align_df.columns.get_loc("TelomereArm")
    qry_label_col = telomeric_align_df.columns.get_loc("QueryLabel")

    
    telomeric_alns = telomeric_align_df.values
    
    # BEGIN NUMPY CODE
    # use numpy array
    
    result_array = np.zeros((len(telomeric_mols), 4))
    
    for i, molecule in enumerate(telomeric_mols):
        # get alignment of specific molecule
        align = telomeric_alns[telomeric_alns[:, qry_id_col] == molecule[0]]
        
        # in our values, our columns for orientation and telomere arm
        # are 2 and 3
         
        orient = align[0,orientation_col]
        telarm = align[0,telomere_col]
        
        # our column for query label is 6
        query_label_idx = int(align[0,qry_label_col])
        
        # nucleotide label position is 
        query_label_nucleotide = molecule[query_label_idx]  
        
        # length of the molecule is the last value
        query_len = molecule[-1]
        
        # get the number of unpaired molecules
        # total number of labels is when you exclude first (ID) and last (LEN)
        # base
        
        # if we get mult
        # LABELPOS-1 | telarm x orientation = 1
        # N-LABELPOS | telarm x orientation -1
        num_labs_query = len(molecule[1:-1])
        # according to formula for telomere length and unpaired labels
        if orient * telarm == 1:
            tlen = query_label_nucleotide
            unpaired_labels = query_label_idx - 1
        if orient * telarm == -1:
            tlen = query_len - query_label_nucleotide
            unpaired_labels = num_labs_query - query_label_idx
        
        result_array[i, :] = [molecule[0], tlen, unpaired_labels,query_len]
    # END NUMPY CODE
    result_df = pd.DataFrame(result_array,columns=["QryContigID","TelomereLen","UnpairedMoleculeLabels","MoleculeLen"])
  
    result_df = pd.merge(telomeric_align_df,result_df,left_on="QryContigID",right_on="QryContigID")
    result_df["Telomere"] = result_df["TelomereArm"].map({1:"left",-1:"right"})
    result_df = result_df.drop("TelomereArm",axis=1)
    
    #reorder cols to be in line with BNGo
    cols = ["RefContigID",
            "QryContigID",
            "MoleculeOrientation",
            "MoleculeConfidence",
            "AlignedLabelPosition",
            "BoundReferencePosition",
            "UnpairedReferenceLabels",
            "UnpairedMoleculeLabels",
            "EndDistance",
            "MoleculeLen",
            "Telomere", # tis is TelomereArm but returned to 'left' 'right'
            "TelomereLen",
            "FullAlignment"
            ]
    return result_df[cols]




# %% FINAL FUNCTIONS - FOR DISPLAYING INFORMATION AND CALCULATING STATISTICS

# GETTING STATISTICS
def get_rough_xmap_stats(xmap_df: pd.DataFrame(), suffix: str = "") -> Tuple[Dict[str, float], Set[int]]:
    """
    Gets rough statistics and returns unique molecules that align to the reference
    Statistics include:
        Length of queried olecule (mean and std)
        Confidence of alignnment (mean and std)
        Number of molecules
        Number of unique molecules
    """
    if len(xmap_df) == 0:
        return dict(), set()



    stats_dict = dict()

    stats_dict[f"qrylen_mean_{suffix}"] = xmap_df["MoleculeLen"].mean()
    stats_dict[f"qrylen_std_{suffix}"] = xmap_df["MoleculeLen"].std()

    stats_dict[f"conf_mean_{suffix}"] = xmap_df["MoleculeConfidence"].mean()
    stats_dict[f"conf_std_{suffix}"] = xmap_df["MoleculeConfidence"].std()

    unique_mols = set(xmap_df["QryContigID"])

    stats_dict[f"molecules_{suffix}"] = len(xmap_df)
    stats_dict[f"unique_mols_{suffix}"] = len(unique_mols) / len(xmap_df)
    
    return stats_dict, unique_mols

def get_xmap_statistics(xmap_df:pd.DataFrame,name:Optional[str]=None) -> Dict[str, float]:
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
    
    xmap_full = xmap_df.loc[xmap_df.FullAlignment==1]
    xmap_part = xmap_df.loc[xmap_df.FullAlignment==0]
    
    #name = os.path.basename(xmap_path).replace(".xmap", "")
    
    stats_full, mols_full = get_rough_xmap_stats(xmap_full, "full_alignment")
    stats_part, mols_part = get_rough_xmap_stats(xmap_part, "partial_alignmet")

    stats = stats_full.copy()
    stats.update(stats_part)
    if mols_part:
        bothmols = len(mols_full.intersection(mols_part))
        jacc = len(mols_full.intersection(mols_part))/len(mols_full.union(mols_part))
        stats["joint_mols"] = bothmols
        stats["jaccard_btw_alns"] = jacc
    for k, v in stats.items():
        logger.info("%s - STATS: '%s': %f", name, k, v)
        
    # return stats
    
def descriptive_telomere_stats(telomere_df:pd.DataFrame,name:Optional[str]=None,
                               groupby:List[str]=["RefContigID","FullAlignment","Telomere"]
                               )->None:
    "Gets descriptive statistics abuot length of telomeres"
    desc = telomere_df.groupby(groupby)["TelomereLen"].describe()
    
    for _, row in desc.iterrows():
        name = dict(zip(["RefContigID","FullAlignment","Telomere"], row.name))
        name = " , ".join([f"{k}:{v}" for k,v in name.items()])+" | "
        display = ' , '.join([f"{index} : {value:.0f}" for index, value in row.items()])
        logger.info("%s - TELOMERE_STATS - %s", name, name + display)

# FILTERING
def reduce_dataset(data: pd.DataFrame, ref_tol: int,
                   mol_tol: int, dis_tol: int) -> pd.DataFrame:
    """
    Reduces our dataset based on the numbers provided.
    For every alignment we have the following three statistics:
        UnpairedReferenceLabels
        UnpairedMoleculeLabels

    Also, additional parameter ('dis_tol') checks 'EndDistance' which
    is the distance between first aligned label and the end of the chromosome.
    In case where our label is too far from the chromosome ends, we can simply
    exclude those molecules from calculation.

    Our alignment is MOLECULE to REFERENCE and this alignment
    doesn't have to be perfect. For example, we can align the penultimate
    label of the CONTIG to last label of REFERENCE.

    These two statistics tell us how many labels are after the last aligned
    pair. Depending on your level of tolerance, this can vary.

    For example, passing 'ref_tol'=1 means that you tolerate one extra label
    after the last paired label on the reference.
    See recommended values in .const.
    """
    input_df = data.copy()

    original_len = len(data)
    len_input_df = original_len

    cols = ["EndDistance", "UnpairedReferenceLabels",
            "UnpairedMoleculeLabels"]
    for name, tol in zip(cols, [dis_tol, ref_tol,  mol_tol]):
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


def calculate_telomere_lengths(xmap_path: str, bnx_path: str, cmap_path: str,
                               partial=True, chunksize:int = 500000,
                               arms:Optional[TelarmType]=["left","right"],
                               pxmap_path:Optional[str]=None,
                               ref_tol:int=REF_TOL,mol_tol:int=MOL_TOL,
                               dis_tol:int=DIS_TOL) -> pd.DataFrame:
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


    logger.info("Getting alignments from XMAP - Full alignments")
    telomeric_alns_df = read_fandom_xmap(xmap_path, cmap_path,how=arms,chunksize=chunksize)
    telomeric_alns_df["FullAlignment"] = 1
    
    # GET PARTIAL ALIGNMENTS
    if pxmap_path is None:
        logger.info("Inferring partial alignment")
    pxmap_path = get_partial_alignment(xmap_path)
    
    if partial and os.path.isfile(pxmap_path):
        logger.info("Getting alignments from XMAP - Partial alignments")
        telomeric_alns_part = read_fandom_xmap(pxmap_path, cmap_path,how=arms,chunksize=chunksize)
        telomeric_alns_part["FullAlignment"] = 0
        telomeric_alns_df = pd.concat([telomeric_alns_df,telomeric_alns_part])
    else:
        logger.info("Skipping partial alignments")
    
    molecule_ids = set(telomeric_alns_df.QryContigID.astype(str).unique())
    
    logger.info("Finding information for %d molecules", len(molecule_ids))
    molecules = read_queried_molecules(bnx_path, molecule_ids)
    
    telomeres = calculate_telomeres(molecules,telomeric_alns_df)
    # filter the molecules
    telomeres = reduce_dataset(telomeres,ref_tol=ref_tol,mol_tol=mol_tol,dis_tol=dis_tol)
    
    # statistics on telomeric_alns_df
    get_xmap_statistics(telomeres,os.path.basename(xmap_path).replace(".xmap", ""))
    
    return telomeres