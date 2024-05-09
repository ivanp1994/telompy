# -*- coding: utf-8 -*-
"""
Created on Thu May  9 12:28:10 2024

@author: ivanp
"""

from telompy.funcs import read_map_file, joinpaths
from telompy.const import CONTIG_PATH, QUERYCMAP_PATH, MASTER_XMAP, MASTER_REFERENCE

_path =  r"C:\#BIONANO_DATA\PROBLEMATIC_TELOMER"


master_xmap = read_map_file(joinpaths(_path, MASTER_XMAP)) #<-the entire contig
master_cmap = read_map_file(joinpaths(_path, MASTER_REFERENCE)) #<-the reference

from telompy.funcs import fish_last_aligned_label

last_aligned = fish_last_aligned_label(master_xmap)
last_aligned = last_aligned.drop_duplicates(["RefContigID", "SiteID", "QryContigID"])
last_aligned.pop("Alignment")

last_label = master_cmap[master_cmap.NumSites==master_cmap.SiteID]
