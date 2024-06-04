# TelOMPy
A Python tool to determine individual (molecule) level of telomere lengths from optical mapping 

[![Licence](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/irahorecka/sgd-rest/main/LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Typing](https://github.com/ivanp1994/telompy/actions/workflows/typehinting.yml/badge.svg)](https://github.com/ivanp1994/telompy/actions/workflows/typehinting.yml)
[![Linting](https://github.com/ivanp1994/telompy/actions/workflows/flaking.yml/badge.svg)](https://github.com/ivanp1994/telompy/actions/workflows/flaking.yml)

# TelOMPy usage for BNGO data
## Input and output parameters
The primary input for `telompy` is a folder with BNGO assembly which is downloaded either from Bionano Access or is the output of Bionano Solve. The structure of said folder (and relevant files) are elucidated in *Structure of BNGO de novo assembly*.

There are two ways to input the BNGO folder:

 - via `--conf`/`-c` parameter which takes a `CSV` (no column names!) file of two columns (second one can be empty) where first column contains paths to BNGO *de novo* assembly folders and the second column contains how the telomere lengths will be saved.

- via `--input`/`-i` and `--name`/`-n` parameters where the first is a white-space delimited list of paths to BNGO *de novo* assembly folders and the is a white-space delimited list of how telomere lengths will be saved.

The output where telomeres will be saved are specified via `--output`/`-o` folder. Telomeres will be saved as `CSV` files in said folder. If there is no `--name`/`-n` parameter or `--conf`/`-c`'s second column is empty (or some elements are empty), then the telomeres will be saved as the base-name of said folder.

Multiprocessing is implemented via `--thread`/`-t` parameter,

Example usage:

`telompy -c conf.csv -o telomeres_bngo -t 8`

and 

`telompy -i BNGO/K55_-_De_novo_pipeline_results BNGO/K56_-_De_novo_pipeline_results -n 5455_Kidney_WD 5456_Kidney_WD -o telomeres_bngo -t 8
`
## Aditional parameters for filtering
There are five additional parameters that are used to filter out the data (for details, see *Theory of operations*):

`--arms`/`-a` - Which ends of chromosomes to calculate telomeres (valid options are 'L' for left and 'R' for right, or 'L' and 'R' for both

`--ref_tol`/`-rt` - Maximum number of unpaired labels on the reference (chromosome) after/before the last aligned pair

`--con_tol`/`-ct` - Maximum number of unpaired labels on the contig (chromosome) after/before the last aligned pair

`--mol_tol`/`-mt` - Maximum number of unpaired labels on the molecule (chromosome) after/before the last aligned pair

`--dis_tol`/`-dt` - Maximum distance between the first/last aligned label on the reference and the end of a chromosome (in nucleotide base pairs).


## Additional parameters for future proofing
These parameters can be changed via CLI, but are also found in `const.py`, and relate to the organizational structure of BNGO folder - for more details, see *Structure of BNGO de novo assembly*.

`-mx`/`--chrom_xmap` - an `XMAP` file of all contigs versus a reference

`-mr`/`--chrom_reference` -a `CMAP` file of the reference (the final, presumably chromosome-level reefrence) 

`-cq`/`--chrom_query` - a `CMAP` file of assembled contigs, serves as an query against `chrom_reference`

`-cf`/`--contig_xmap_format` - an `XMAP` file of queried molecules against a *de novo* assembly of queried molecules

`-qc`/`--contig_query_format` - a `CMAP` file of molecules that serves as a query, and first step of assembly


## Installation and usage

The only thing required is pandas, as outlined in `requirements.txt`. 
Additionally, a singularity/apptainer definition file is given in `telompy.def`. 

It can be used as a Python API (the main function being `calculate_telomere_lengths`) which takes all the parameters the CLI tool takes,
and it can be used from the CLI.

# Theory of operation

We aim to provide a relatively sensible way to determine the telomere lengths. 

Telomeres are repetitive nucleotide sequences associated with specialized proteins at the ends of linear chromosomes commonly found in eukaryotes.
Optical mapping, as implemented by Bionano Genomics (BNGO) is a technique that involves marking high molecular weight (HMW) DNA on specific motifs, and then 
assembling a complete map of an organism from the patterns of said marking.

First off, we actually cannot determine *absolute* telomere length. Telomeres are repetitive sequences, and specific motifs recognized by BNGO enzymes cannot mark a telomere.
What we can do, is determine the end point after (for the right end of  thechromosome) or before (for the left end of the chromosome) the last or the first label on the reference.
Visualized below:

![drawing-1](https://github.com/ivanp1994/telompy/assets/84333373/d1ae2dd4-6230-495f-a136-78d67b475a66)


The green rectangles represent a reference, whilst blue rectangles represent a single molecule. Optical mapping works by aligning labels (thin blue lines) on a molecule to 
labels on the reference. The red squares are segments of molecule before the first label that aligns to the first label on the reference and after the last label that aligns to the last label on the reference.
We can call this *relative* telomere length.
The procedure is simple - we find the bound label (last label for the right, and first label for the left) on the reference, we find its pair on the individual molecule,
and then (depending on the orientation of the assembly), the *relative* telomere length is the part of the individual molecule *after* or *before* said pair.

The procedure is made more complicated by the way BNGO assembles and maps. In the first round, the individual molecules are assembled into contigs.
In the second round, the assembled contigs are aligned to the reference (or de novo assembled if no reference is given). This gives two distinct views:

  1) Molecules to contigs
  2) Contigs to reference

So in order to find molecules mapped to reference, we must *overlay* these two views to get "molecules-to-reference".
One problem that can occur is related to high coverage (and therefore a huge amount of individual molecules needed),
the reason why high coverage of optical mapping is necessary (large amounts of false positives) and possible mutations in individual molecules
that can plop a label where there should be none.

For example:

![15_101](https://github.com/ivanp1994/telompy/assets/84333373/1e162922-3b25-4f14-b4e5-39aa798e7ec6)


In the above image, there is a label on contig (marked yellow) that is not paired to the label of reference - the penultimate label on the contig maps
to the last label on the reference. What happened here is that out of a bunch of molecules that assembled said contig, a sufficient number of molecules had
one extra label - a mutation or a false positive (an enzyme labeled a place it shouldn't have) happened on - so our contig got one extra label that does not fit in the reference.

Similar problems can occur for a reference - if the **last/first aligned** label on the reference is not the **last/first** label on the reference, we cannot call this a *relative* telomere, let alone a telomere.
To control for this, we implement a few key filters elucidated in *Additional parameters for filtering*. 

Recapped, those are related to the maximum number of labels on reference/contig/molecule before/after the first/last **aligned** pair of molecule-contig-reference.
Note the bolded part - if you have e.g. 3 extra labels on the contig that do not align to the reference, and have 15 extra labels on the molecule that do not align to the reference, but 3 of those actually align to the contig,
we are still counting 15 extra labels.

One problem that can also occur:

![1_lchrom](https://github.com/ivanp1994/telompy/assets/84333373/3779f077-4f82-4561-a5e2-6802ed011cef)


In here, our first label is found good 3 million bases after the chromosome start. The first strech of the reference is unlabeled.
We cannot reasonably call this a telomere - and to that extend we define the maximum distance between the first/last **aligned** label on the reference and the end of the chromosome it's on.
So an additional parameter `dis_tol` is here to control that - the exact value of it depends on the annotation and the assembly of your model organism.
For mouse, as visualized above, left telomeres are annotated very poorly, and we cannot reasonably calculate them, even as a *relative* value. 


# TelOMpy output structure

The file is saved as a `CSV` file with the following columns:

`RefContigID` - the ID of reference contig (the chromosome)

`MoleculeID` - the ID of a molecule

`QryContigID` - the ID of an assembled contig

`MoleculeOrientation` - the orientation of Molecule-Contig assembly (1 for + and -1 for -)

`ContigOrientation` - the orientation of Contig-Reference assembly (1 for + and -1 for -)

`MoleculeConfidence` - the confidence score of Molecule-Contig assembly

`AlignedLabelPosition` - the position (in nucleotide basepairs) of the first/last aligned label on the reference

`BoundReferencePosition` - the position (in nucleotide basepairs) of the first/last label on the reference

`UnpairedReferenceLabels` - the number of labels on the reference before/after the last aligned label

`UnpairedContigLabels` - the number of labels on the contig before/after the last aligned label

`UnpairedMoleculeLabels` - the number of labels on the molecule before/after the last aligned label

`EndDistance` - the distance (in nucleotide basepairs) between the start/end of the chromosome and first/last aligned label

`MoleculeLen` - the length of the molecule

`Telomere` - left or right telomere

`TelomereLen` - the length of the telomere

- - -
Some notes:

"AlignedLabelPosition" and "BoundReferencePosition" will be equal if "UnpariedReferenceLabels" is 0.

"MoleculeConfidence" is defined as "Negative Log10 of p-value of alignment" and details can be found in Bionano's Theory of Operations - Structural Variant Calling pdf document.
I'll C/P what I've learned from BNGO's support:
>There is a general description on page 17, and more in the appendix sections for different variant types. Please also find a more detailed description for the p-values for insertions and deletions confidence scores below and in the attached pdf (publication on which the modeling is based):
> The non-normalized likelihood ratio for the insertion or deletion region is based on Gaussian error likelihood and the ratio of label interval vs Gaussian SD
> The confidence is based on the average PPV (fraction of calls that are correct) for a large number of simulated insertions (or deletions) in the same size bin (eg 500-1000bp, 1kb-2kb etc) and same Pvalue or Likelihood ratio bin (whichever is less stringent, typically that is the Gaussian likelihood ratio, unless the insertion or deletion is near the end of the alignment, and the p-Value bins are 0.1+, 0.01-0.1, 0.001-)
> **A confidence score of 0 could also indicate an area of high (or low) complexity, for which alignments were difficult to determine. It does not necessarily mean a variant should not be considered at all.**  

# TelOMpy logging structure

Other than finished `CSV` file, `telompy` also writes log to output. Output is in the form of `[LEVEL] message`.
One input cointains `LABEL_STATS` flag which gives out the number of unpaired labels on the reference and contig level as well as the distance of the end of the chromosome to the relevant label.
For example:

```
[INFO] 03-1_-_De_novo_pipeline_results - LABEL_STATS - RefContigID : 2 , QryContigID : 22 , UnpairedReferenceLabels : 0 , UnpairedContigLabels : 0 , EndDistance : 3060417
```

On chromosome number 2 for contig number 22, we have 0 unpaired references on both the contig and the chromosome, but our distance between the end of the chromosome and the aligned label (`EndDistance`) is roughly 300kbp.

Additionally, the log includes detailes about filtering in the form of number of molecules excluded for that particular filter:

```
[INFO] Excluded EndDistance > 500000 from input_df - went from 3266 to 57 molecules (98.25 percent reduction)
[INFO] Excluded UnpairedReferenceLabels > 0 from input_df - went from 57 to 0 molecules (100.00 percent reduction)
[INFO] Went from 3266 to 0 molecules (100.00 percent reduction)
```

or:

```
[INFO] Excluded EndDistance > 500000 from input_df - went from 2520 to 2444 molecules (3.02 percent reduction)
[INFO] Excluded UnpairedReferenceLabels > 0 from input_df - went from 2444 to 2233 molecules (8.63 percent reduction)
[INFO] Excluded UnpairedContigLabels > 8 from input_df - went from 2233 to 2233 molecules (0.00 percent reduction)
[INFO] Excluded UnpairedMoleculeLabels > 200 from input_df - went from 2233 to 2233 molecules (0.00 percent reduction)
[INFO] Went from 2520 to 2233 molecules (11.39 percent reduction)
```

The next stats that are thrown at the stdout are:
1. Correlation between molecule lengths and telomere lengths
2. Z-testing between molecule lengths of telomeric molecules and expected molecule length

Number one relates to the fact that the bigger the molecule, the bigger expected telomere length. We are using `pearsonr` with no correction for multiple testing as every element is calculated independently.
Number two checks if the telomeric molecules differ significantly in length from the average expected molecule length. Expected molecule length is calculated as average molecule length of all molecules that 
enter aligning process and this information can be found in `exp_informaticsReportSimple.json` file.

Example:
```
[INFO] 03-1_-_De_novo_pipeline_results - LABEL_STATS - MoleculeLen/TelomereLen correlation 0.4464 with pvalue 0.00000000 
[INFO] 03-1_-_De_novo_pipeline_results - LABEL_STATS - Difference from expected molecule length -  Cohen D is -0.2380 with pvalue 0.00000000 
```


# Structure of BNGO *de novo* assembly

The (relevant) structure of an un-zipped BNGO assembly is as follows:

```
output/
└── exp_informaticsReportSimple.json
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
(If one wants to use alternative optical map assembly tools, such as [FaNDOM](https://github.com/jluebeck/FaNDOM) the (processed) molecules are given in `autoNoise1_rescaled.bnx` and the reference is given in `exp_refineFinal1_merged_q.cmap`)

The BNGO assembly is done in two parts:
First a *de novo* assembly of contigs from molecules where molecules are assembled into N contigs of specific IDs. The contigs are found in the folder `refine1_ExperimentLabel` in the three forms:

 - EXP_REFINEFINAL1_contig{x}.xmap - which is an `XMAP` file of a segment of molecules versus *de novo* contigs
 - EXP_REFINEFINAL1_contig{x}_r.cmap - which is a `CMAP` file of a reference contig
 - EXP_REFINEFINAL1_contig{x}_q.cmap - which is a `CMAP` file of the raw molecules

The `{x}` is a specific contig ID.

The second step is to align contigs to the reference. The files are found in a folder above `refine1_ExperimentLabel` and consist of three important files:
- exp_refineFinal1_merged.xmap - which is an `XMAP` file of all contigs versus a reference
- exp_refineFinal1_merged_r.cmap - which is a `CMAP` file of our reference
- exp_refineFinal1_merged_q.cmap - which is a `CMAP` file of the query contigs.

Since optical mapping by BNGO can easily achieve chromosome-level coverage, we will refer to the first step as "molecule-to-contig" and the second step as "contig-to-chromosome" alignment.
All locations of relevant files are specified in `const.py` for (relative) future-proofing should BNGO deem a change in the structure is necessary.

