# -*- coding: utf-8 -*-
"""
Created on Wed May  8 13:51:38 2024

@author: ivanp
"""
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

path =  r"C:\#BIONANO_DATA\PROBLEMATIC_TELOMER"
path = r"E:\isabella\telomeres\PROBLEMATIC_TELOMER"


telomeres_no_gap = calculate_telomere_lengths(path,gap_size=None)

telomeres_gap = calculate_telomere_lengths(path,gap_size= 100000)

#%%
tel_nogap = pd.concat(telomeres_no_gap)
tel_gap = pd.concat(telomeres_gap)

#%%
minstat_nogap = tel_nogap.groupby("RefContigID")["TelomereLen_corr"].min()
minstat_gap = tel_gap.groupby("RefContigID")["TelomereLen_corr"].min()


#%% the above proves that when I include gap size, I get negative telomeres,
# but negative telomeres are in chrom
# 21 20 16 12 10 7 6 5 4 3 2 1

# basically everywhere!!!!

#further exploration is required!