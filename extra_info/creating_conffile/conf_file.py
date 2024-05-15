# -*- coding: utf-8 -*-
"""
Created on Wed May 15 13:14:03 2024

@author: ivanp
"""
import pandas as pd
import logging
logger = logging.getLogger()
import numpy as np
#open paths
with open("zip_files.txt") as _f:
    fajls = [f.strip().replace(".zip","") for f in _f.readlines()]



excel = pd.read_excel("ProjectOverview.xlsx",skiprows=1)
excel = excel.loc[excel.Experiment=="WD"]

excel["ID"] = excel.MouseID.astype(str) + "_" + excel.Tissue + "_" + excel.Treatment
mapper = dict(zip(excel["LabId"],excel["ID"]))


confa = pd.DataFrame(fajls,columns=["FOLDER_NAME"])
confa["LabId"] = confa["FOLDER_NAME"].str.split('_').str[0] 

confa["name"] = confa.LabId.map(mapper)
confa.pop("LabId")

confa.FOLDER_NAME = "/beegfs/home/ipokrova/BNGO/" + confa.FOLDER_NAME

confa.to_csv("conf.csv",index=False,header=False)



def target_from_config(args):
    "constructs target from a given config file"
    if args["conf"] is None:
        return None
    try:
        conf = pd.read_csv(args["conf"], header=None)
    except FileNotFoundError:
        logger.error("No configuration file found at '%s'", args["conf"])
        return None

    if conf.shape[1] < 2:
        logger.error("Two columns must be in configuration file")
        return None
    # E1101 (no-member) : Instance of 'TextFileReader' has no 'replace'
    # this is a false positive
    conf = conf.replace({np.nan: None})  # pylint: disable=E1101
    return list(zip(conf[0], conf[1]))

fila = target_from_config(dict(conf="conf.csv"))
