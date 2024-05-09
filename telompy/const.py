# -*- coding: utf-8 -*-
"""
Created on Wed May  8 10:19:03 2024

@author: ivanp
"""
import logging

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)

LOGGER = logging.getLogger("telompy")
MASTER_XMAP = "output/contigs/annotation/exp_refineFinal1_merged.xmap"
MASTER_QUERY = "output/contigs/annotation/exp_refineFinal1_merged_q.cmap"
MASTER_REFERENCE = "output/contigs/annotation/exp_refineFinal1_merged_r.cmap"
CONTIG_PATH = "output/contigs/annotation/refine1_ExperimentLabel/EXP_REFINEFINAL1_contig{x}.xmap"
REFCMAP_PATH = "output/contigs/annotation/refine1_ExperimentLabel/EXP_REFINEFINAL1_contig{x}_r.cmap"
QUERYCMAP_PATH = "output/contigs/annotation/refine1_ExperimentLabel/EXP_REFINEFINAL1_contig{x}_q.cmap"
