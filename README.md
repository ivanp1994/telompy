This branch performs telomere length calculations equivalent to those found in `main` (or `left`) branch but for alignmeents made by [FaNDOM](https://github.com/jluebeck/FaNDOM).
FaNDOM differs from BNGO's assembly in that the output is an `XMAP` file of molecules-to-reference - there is no intermediary steps of assembling contigs and therefore no contig level to account for.
This is good as it eases the calculation, and it's bad as the resulting `XMAP` file is large, and the input molecules are also large so memory issues can arise here.

As a reminder the structure of *de novo* BNGO assembly is:
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
So the input can be just the location of the folder. However, input in this case must be a pointer to a `BNX` file (query), a refrence `CMAP` file (reference), and an `XMAP` file that's the result of
FaNDOM assembly procedure.
These inputs can be given via `-bnx`, `-ref`, `-xmap` parameters, or can be given via `-c`/`--conf` file. In the latter case, the `-c`/`--conf` file must be a CSV (no header) file of a minimum of 
two columns (`xmap` and `bnx`) with third being `ref` (or if the `ref` is the same between all mappings, it can be provided via `-ref` parameter in conjunction with `-c`/`--conf` file), and fourth
being how telomeres will be saved.

Filtering options are the same as the main branch, with the exception of `con_tol` as there are no contigs.

The resulting output is identical with as the output of the main branch, with an additional column `FullAlignment`.
FaNDOM aligns molecules in two steps - a "full" alignment and a "partial" alignment (details on what is considered partial and what is considered full alignment can be found in their paper).
Usually this means two `XMAP` files with one being suffixed with `_partial.xmap`. `FullAlignment` contains information if the assembly was from full (1) or partial (0) alignment.

As of yet, this branch is untested and there are problems with type hinting.
