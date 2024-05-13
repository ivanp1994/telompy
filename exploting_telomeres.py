# -*- coding: utf-8 -*-
"""
Created on Mon May 13 13:25:06 2024

@author: ivanp
"""
import seaborn as sns
import os
from telompy import calculate_telomere_lengths


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

PATH =  r"C:\#BIONANO_DATA\PROBLEMATIC_TELOMER"
#path = r"E:\isabella\telomeres\PROBLEMATIC_TELOMER"



def read_telomeres(path=PATH,gap_size=100000):
    if os.path.isfile("telomere_lengths/telomeres.csv"):
        telomeres_gap = pd.read_csv("telomere_lengths/telomeres.csv").drop_duplicates("MoleculeID")
        return telomeres_gap
    telomeres_gap = calculate_telomere_lengths(path=path,gap_size=gap_size) #it returns a list -.-
    pd.concat(telomeres_gap,ignore_index=True).to_csv("telomere_lengths/telomeres.csv",index=False)

    return read_telomeres(path=path,gap_size=gap_size)


telomeres = read_telomeres()
telomeres = telomeres.drop_duplicates("MoleculeID")
data = telomeres.copy()
#%%
above_zero =  data.groupby('RefContigID')['TelomereLen_corr'].apply(lambda x: (x > 0).sum())
below_zero =  data.groupby('RefContigID')['TelomereLen_corr'].apply(lambda x: (x < 0).sum())

