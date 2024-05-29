# TeOMPy usage for BNGO data
## Input
The primary input for `telompy` is a folder with BNGO assembly which is downloaded either from Bionano Access or is the output of Bionano Solve. The structure of said folder (and relevant files) are elucidated in *Structure of BNGO de novo assembly*.

There are two ways to input the folder:
 - 1 via `--conf`/`-c` parameter:
The parameter takes a `CSV` (no column names!) file of two columns (second one can be empty) where first column contains paths to BNGO *de novo* assembly folders and the second column contains how the telomere lengths will be saved.
- 2 via `--input`/`-i` and `--name`/`-n` parameters where the first is a white-space delimited list of paths to BNGO *de novo* assembly folders and the is a white-space delimited list of how telomere lengths will be saved.

Example usage (via apptainer):

`telompy -c conf.csv -o telomeres_bngo -t 8`

and 

`telompy -i BNGO/K55_-_De_novo_pipeline_results BNGO/K56_-_De_novo_pipeline_results -n 5455_Kidney_WD 5456_Kidney_WD -o telomeres_bngo -t 8
`


# Structure of BNGO *de novo* assembly

The (relevant) structure of an un-zipped BNGO assembly is as follows:

```
output/
└── contigs/
    └── auto_noise/
        └──autoNoise1_rescaled.bnx
    └── annotation/
        ├── exp_refineFinal1_merged.xmap
        ├── exp_refineFinal1_merged_q.cmap
        ├── exp_refineFinal1_merged_r.cmap
        └── refine1_ExperimentLabel/
            ├── EXP_REFINEFINAL1_contig{x}.xmap
            ├── EXP_REFINEFINAL1_contig{x}_r.cmap
            └── EXP_REFINEFINAL1_contig{x}_q.cmap
```
`XMAP` (*cross*-map) files represent an alignment file (equivalent to a `BAM` or `SAM` file) between a query and a reference.
Reference must be in `CMAP` (*consensus*-map) format, while the query can either be in `BNX` or `CMAP` format.
The reference query in `BNX` format is found as `autoNoise1_rescaled.bnx`, but we utilize reference queries in the `CMAP` format found either as `exp_refineFinal1_merged_q.cmap` or ˛`EXP_REFINEFINAL1_contig{x}_q.cmap`.

The BNGO assembly is done in two parts:
First a *de novo* assembly of contigs from molecules where molecules are assembled into N contigs of specific IDs. The contigs are found in the folder `refine1_ExperimentLabel` in the three forms:

 - EXP_REFINEFINAL1_contig{x}.xmap - which is an `XMAP` file of a segment of molecules versus *de novo* contigs
 - EXP_REFINEFINAL1_contig{x}_r.cmap - which is a `CMAP` file of a reference contig
 - EXP_REFINEFINAL1_contig{x}_q.cmap - which is a `CMAP` file of the raw molecules

The second step is to align contigs to the reference. The files are found in a folder above `refine1_ExperimentLabel` and consist of three important files:
- exp_refineFinal1_merged.xmap - which is an `XMAP` file of all contigs versus a reference
- exp_refineFinal1_merged_r.cmap - which is a `CMAP` file of our reference
- exp_refineFinal1_merged_q.cmap - which is a `CMAP` file of the query contigs.

Since optical mapping by BNGO can easily achieve chromosome-level coverage, we will refer to the first step as "molecule-to-contig" and the second step as "contig-to-chromosome" alignment.
All locations of relevant files are specified in `const.py` for (relative) future-proofing should BNGO deem a change in the structure is necessary.

# Steps to calculate (relative) telomere length
