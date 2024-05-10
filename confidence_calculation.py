# -*- coding: utf-8 -*-
"""
Created on Thu May  9 12:28:10 2024
TelomerePosition!

we can rip telomere positions from Table Tool
https://genome.ucsc.edu/cgi-bin/hgTables


This information can be found in the "gap" database table. Use the Table Browser to extract it. To do this, select your assembly and the gap table, then click the " filter Create" button. Set the "type" field to centromere telomere (separated by a space). For help using the Table Browser, visit the User's Guide.


What we find here is following:
    
    For right telomeres (start is not 0, using start)
    
    ContigLength - TelomerePosition is equal to gap size
    
    In right telomeres, gap size is 100000 (arbitrary), meaning that
    our length of contig is equal to TelomereStart + Gap Size
    
    TELOMERE:   ------------|-----------]
    REFERENCE:  ----x-------------------]
    
    in the above, x is the last label on reference
    
    We can estimate distance of Telomere to Label pos
    
    TELOMERE:   ------------|-----------]
    REFERENCE:   ----x-------------------]
    DISTANCE:  ------++++++++-------------
    
    SINCE OUR BASICALLY WORKS LIKE THIS
    
    REFERENCE: ---------|------------]
    MOLECULE:  ---------|----|-----|-----]
    TELOMERE:           |----------------|
    
    we would theoretically need to subtract! the distance
    from our telomere length to get better calculation
    
    
More details

x = Labels
] = End of contig
    

REFERENCE:  ---------x-----]
    
We know that the gap size is 10 000, lets call this GS
Our position of telomere start is equal to:
    RefLen - GS
    
From this, we have two options:

1. TelomereStart (TS) is behind the label:
REFERENCE: ---------|--x-----]

2. TelomereStart(TS) is in front of the label:
REFERENCE: -----------x---|--]

In first case, TS (|) - LL (x) is going to be
negative     
so we can theoretically calculate Telomere_corr via

Telomere_corr = TelomereLen - Offset




    
    
    
    THEORETICALLY
    
TelomerePos = RefLen - 100000 #gap size
TS = TelomereStart
RL = ReferenceLength
LL = LastLabel
GS = GapSize
O = Offset

TL = TelomereLength
TLC = TelomereLengthCorrected

TS = RL - GS 
O = LL - TS

O = LL - (RL - GS)
O = LL - RL + GS


TLC = TL - O
    


ContigLen == RefLen

@author: ivanp
"""
#%%
import pandas as pd
from telompy.funcs import read_map_file, joinpaths
from telompy.const import CONTIG_PATH, QUERYCMAP_PATH, MASTER_XMAP, MASTER_REFERENCE, MASTER_QUERY
from telompy.funcs import fish_last_aligned_label,fish_last_label
from telompy.funcs import read_molecules_cmap
from telompy.funcs import get_first_query_last_reference
from telompy.funcs import read_contig_xmap
from telompy.funcs import get_number_of_unpaired_contig_labels
from telompy.funcs import calculate_telomere

from typing import Literal,Union

path =  r"C:\#BIONANO_DATA\PROBLEMATIC_TELOMER"
path = r"E:\isabella\telomeres\PROBLEMATIC_TELOMER"

#master_xmap = read_map_file(joinpaths(path, MASTER_XMAP)) #<-the entire contig
#master_cmap = read_map_file(joinpaths(path, MASTER_REFERENCE)) #<-the reference
#master_query = joinpaths(path,MASTER_QUERY,)



#read in master reference
master_cmap = read_map_file(joinpaths(path, MASTER_REFERENCE))
length = master_cmap.drop_duplicates(["CMapId","ContigLength"])[["CMapId","ContigLength"]]






#%% now we have a file called LastPosition in MASTER XMAP
master_xmap = fish_last_label(joinpaths(path))
master_xmap.pop("Alignment")

#%% load in telomere positions 
_telomerepos = pd.read_csv(r"C:\Users\Ivan\Desktop\IRB_Folder\TELOMERE\telomerepos.txt",
            sep="\t")

_telomerepos["chrom"]=_telomerepos["chrom"].str.replace("chr","").replace({'X': 20, 'Y': 21}).astype(int)
teloms_right = _telomerepos[_telomerepos["chromStart"]!=0][["chrom","chromStart"]]
teloms_right.columns = ["RefContigID","TelomerePosition"]


master_xmap = pd.merge(master_xmap,teloms_right,
                       left_on="RefContigID",right_on="RefContigID",
                       how="left")

master_xmap = pd.merge(master_xmap,length,
                       right_on="CMapId",
                       left_on="RefContigID",
                       how="left")

#
master_xmap["RefLen/ContigLength"] = master_xmap["RefLen"] / master_xmap["ContigLength"]



## formula


#%%
master_xmap["Length-TelomerePos"] = master_xmap["ContigLength"] - master_xmap["TelomerePosition"]
master_xmap["TelomerePos-LastAlignedLabel"] = master_xmap["TelomerePosition"] - master_xmap["AlignedLabelPosition"]
master_xmap["TelomerePos-LastLabel"] = master_xmap["TelomerePosition"] - master_xmap["RefEndPos"]

master_xmap["TeloPosCalc"] = master_xmap["RefLen"] -  100000
master_xmap["OffSet"] = master_xmap["TeloPosCalc"] - master_xmap["AlignedLabelPosition"]

gap_size = 10000
master_xmap["OffSet2"] = master_xmap["AlignedLabelPosition"] - master_xmap["RefEndPos"] + gap_size

gap_size = 0
master_xmap["OffSet3"] = master_xmap["AlignedLabelPosition"] - master_xmap["RefEndPos"] + gap_size

#row= master_xmap.iloc[7] # 4th 2611 map
#telamer = calculate_telomere(path=path,row=row)
#%%

#%%
#import seaborn as sns

#sns.relplot(x=telamer["TelomereLen_corr"],y=telamer["UnpairedMoleculeLabels"])


"""



row= master_xmap.iloc[7] # 4th 2611 map
#row = master_xmap.iloc[39] # 20th 1430 map #44 unpaired labels

contig_format = CONTIG_PATH
querycmap_format = QUERYCMAP_PATH

master_contig = row["QryContigID"]
contig_alignment = row["Alignment"]

contig_path = joinpaths(path, contig_format.format(x=master_contig))
molecules_path = joinpaths(path, querycmap_format.format(x=master_contig))

aligned_label = get_first_query_last_reference(contig_alignment)
# reads only those alignments that contain query label
# on the reference
contig_aligned = read_contig_xmap(contig_path, aligned_label) #dobro je

molecules = read_molecules_cmap(molecules_path, contig_aligned.QryContigID.unique())
from telompy.funcs import extract_reference_query_pair


paired_query = extract_reference_query_pair(contig_aligned,aligned_label)
contig_aligned["RefContigMolPair"] = paired_query.astype(str).apply(lambda x: f"({aligned_label},{x})")






contig_aligned["UnpairedMoleculeLabels"] = contig_aligned.apply(_gnoucl,axis=1,
                                                                molecules=molecules)
"""


"""
subrow = contig_aligned.iloc[0] #Orientation is -
subrow = contig_aligned.iloc[1] #Orientation is +

unpara = get_number_of_unpaired_contig_labels(path=molecules,
                                              alignment=subrow["RefContigMolPair"],
                                              qrycontigid=subrow["QryContigID"],
                                              orientation=subrow["Orientation"])


"""


#%% unpara = 6 qryid = 7975902
#unapram = molecules.loc[molecules.CMapId==7975902] # last SITE ID 31
#pair of alignment
# 38, 25
# 31 - 25

#%% CALCULATING TELOMERES




#%% geting unpaired contig labels
"""


last_label_qry = get_first_query_last_reference(row["Alignment"])
"""
def get_number_of_unpaired_contig_labels(path:str,alignment:str, qrycontigid:int,
                                     orientation:Literal["-","+"],                                     
                                     contig_query:str=MASTER_QUERY)->int:
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
    if orientation not in ["-","+"]:
        raise ValueError("Invalid orientation - must be either '-' or '+'")
        
    last_label_qry = get_first_query_last_reference(alignment)
    
    #if orientation is negative
    #then queries are descending
    #so it always ends with 1
    if orientation == "-":
        return last_label_qry - 1
    
    #if it's positive, we need to read molecules cmap
    _cmap_contig_ref = read_molecules_cmap(joinpaths(path,contig_query),[qrycontigid])
    #it's the penultimate SiteID
    #last site ID is contig length
    last_label = _cmap_contig_ref["SiteID"].max()-1
    return last_label - last_label_qry
    

"""
numba = get_number_of_unpaired_contig_labels(path = path,
                                         alignment = row["Alignment"],
                                         qrycontigid = row["QryContigID"],
                                         orientation = row["Orientation"],
                                         contig_query = MASTER_QUERY
                                         )


"""

"""
# this is the number of labels on contig that are not paired
# with reference - after the last pair is aligned
if row["Orientation"]=="-":
    unpaired_contig_labels = last_label_qry-1 
else:
    #reading info
    #last label is the SiteID-1 because
    #site ID last is contig length
    unpaired_contig_labels = read_molecules_cmap(joinpaths(_path,MASTER_QUERY,),
                                      [row["QryContigID"]])["SiteID"
                                                            ].max()-1
    unpaired_contig_labels = unpaired_contig_labels - last_label_qry                                  

"""

#%% FISHING NUMBERS OF UNPAIRED LABELS ON REFERENCE - DONE!!!
#<-last aligned
#<-reference_cmap

"""
from telompy.funcs import fish_last_aligned_label

last_aligned = fish_last_aligned_label(master_xmap)

reference_cmap = master_cmap

aligned_positions = last_aligned[["RefContigID", "SiteID", "QryContigID"]].drop_duplicates()

aligned_positions = pd.merge(left=aligned_positions,
                             right=reference_cmap[["CMapId", "SiteID", "Position"]],
                             left_on=["RefContigID", "SiteID"],
                             right_on=["CMapId", "SiteID"],
                             )[["RefContigID", "QryContigID", "SiteID", "Position"]]

# get positions of last labels
# last label on reference
last_labels = reference_cmap[reference_cmap.SiteID == reference_cmap.NumSites][["CMapId","SiteID","Position"]]
last_labels.columns = ["RefContigID","LastID" ,"LastPosition"]

# merge the two
aligned_positions = pd.merge(left=aligned_positions,
                             right=last_labels,
                             left_on="RefContigID",
                             right_on="RefContigID",
                             how="left")

aligned_positions["Offset"] = aligned_positions["LastPosition"] - aligned_positions["Position"]
aligned_positions["Offset_Label"] = aligned_positions["LastID"] - aligned_positions["SiteID"]
"""