# TeOMPy usage for BNGO data
## Input and output
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
These parameters can be changed via CLI, but are also found in `const.py`, and relate to the organizational structure of BNGO folder.
TODO: dont really want to bother with this


# Theory of operation

We aim to provide a relatively sensible way to determine the telomere lengths. 

Telomeres are repetitive nucleotide sequences associated with specialized proteins at the ends of linear chromosomes commonly found in eukaryotes.
Optical mapping, as implemented by Bionano Genomics (BNGO) is a technique that involves marking high molecular weight (HMW) DNA on specific motifs, and then 
assembling a complete map of an organism from the patterns of said marking.

First off, we actually cannot determine *absolute* telomere length. Telomeres are repetitive sequences, and specific motifs recognized by BNGO enzymes cannot mark a telomere.
What we can do, is determine the end point after (for the right end of  thechromosome) or before (for the left end of the chromosome) the last or the first label on the reference.
Visualized below:
**ADD IMAGE**

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
**ADD IMAGE**

In the above image, there is a label on contig (marked yellow) that is not paired to the label of reference - the penultimate label on the contig maps
to the last label on the reference. What happened here is that out of a bunch of molecules that assembled said contig, a sufficient number of molecules had
one extra label - a mutation or a false positive (an enzyme labeled a place it shouldn't have) happened on - so our contig got one extra label that does not fit in the reference.

Similar problems can occur for a reference - if the **last/first aligned** label on the reference is not the **last/first** label on the reference, we cannot call this a *relative* telomere, let alone a telomere.
To control for this, we implement a few key filters elucidated in *Additional parameters for filtering*. 

Recapped, those are related to the maximum number of labels on reference/contig/molecule before/after the first/last **aligned** pair of molecule-contig-reference.
Note the bolded part - if you have e.g. 3 extra labels on the contig that do not align to the reference, and have 15 extra labels on the molecule that do not align to the reference, but 3 of those actually align to the contig,
we are still counting 15 extra labels.

One problem that can also occur:

**ADD IMAGE**

In here, our first label is found good 3 million bases after the chromosome start. The first strech of the reference is unlabeled.
We cannot reasonably call this a telomere - and to that extend we define the maximum distance between the first/last **aligned** label on the reference and the end of the chromosome it's on.
So an additional parameter `dis_tol` is set. 





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
(If one wants to use alternative optical map assembly tools, such as [FaNDOM](https://github.com/jluebeck/FaNDOM) the (processed) molecules are given in `autoNoise1_rescaled.bnx` and the reference is given in `exp_refineFinal1_merged_q.cmap`)

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

