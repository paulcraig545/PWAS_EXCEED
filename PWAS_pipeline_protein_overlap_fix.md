# PWAS Pipeline Protein List Compatibility Fix

## Overview

This update addresses an error encountered when running `PWAS_pipeline.r` using the EXCEED protein file and the supplied LASSO protein list.

The pipeline reads in:

- a protein matrix/dataframe from `prot_filepath`
- a list of selected LASSO proteins from `probe_filepath`

The original script assumed that every protein in the LASSO protein list was present as an exact column name in the protein file. However, this was not true for the EXCEED protein dataset.

## Error encountered

The job failed with:

```r
Error in [.data.frame(ms_prot, , prot_list) : 
  undefined columns selected
```

This occurred during the protein transformation step, where the script attempted to subset the protein dataframe using the full `prot_list`.

## Cause of the problem

The original script created `ms_prot` using only proteins that overlapped between the protein file and the LASSO protein list:

```r
ms_prot <- proteins %>%
    rename(id = all_of(id_col)) |>
    select(id, intersect(names(proteins), prot_list))
```

However, the script later still used the original full `prot_list`:

```r
prot_std <- mapply(
  transform,
  ms_prot[, prot_list]
)
```

This caused an error because some proteins in `prot_list` were not present in `ms_prot`.

In this specific case, the LASSO protein list contained 19 proteins. Of these, 17 were exact matches in the EXCEED protein file. One protein was only present as part of a compound protein column, and one protein was not present as an exact match.

## Fix applied

The fix explicitly identifies which proteins from the LASSO list are present in the protein file and updates `prot_list` so that all downstream steps use only available proteins.

```r
proteins <- readRDS(prot_filepath)
prot_list <- readRDS(probe_filepath)

# Keep only proteins that exist in the protein file
available_prots <- intersect(colnames(proteins), prot_list)
missing_prots <- setdiff(prot_list, colnames(proteins))

print(paste0("Proteins requested: ", length(prot_list)))
print(paste0("Proteins found in protein file: ", length(available_prots)))
print(paste0("Proteins missing from protein file: ", length(missing_prots)))

if (length(missing_prots) > 0) {
  print("Missing proteins:")
  print(missing_prots)
}

# Update prot_list so downstream code only uses available proteins
prot_list <- available_prots

ms_prot <- proteins %>%
    rename(id = all_of(id_col)) |>
    select(id, all_of(prot_list))
```

## What the fix does

This fix:

1. Reads the protein file and selected LASSO protein list.
2. Identifies proteins that are present in both datasets.
3. Reports how many proteins were requested.
4. Reports how many proteins were found in the protein file.
5. Prints any missing proteins.
6. Updates `prot_list` so downstream analysis only uses valid protein columns.
7. Prevents the `undefined columns selected` error.

## Why this is appropriate

The LASSO weights were trained on a specific set of proteins, but the replication protein file may not contain every protein as an exact column match. Some proteins may be absent, renamed, or represented as compound UniProt IDs such as:

```r
P0DOY2;P0DOY3
```

Using only exact matches avoids forcing an uncertain mapping between a single protein ID and a compound protein column.

## Important methodological note

This fix uses exact protein ID matching only.

For example, if the LASSO list contains:

```r
P0DOY3
```

but the protein file contains:

```r
P0DOY2;P0DOY3
```

this will be treated as missing because it is not an exact column-name match.

This is conservative and avoids introducing assumptions about whether a compound protein measurement should be used as a proxy for an individual protein.

## Expected log output

After applying the fix, the log should include lines like:

```r
Proteins requested: 19
Proteins found in protein file: 17
Proteins missing from protein file: 2
Missing proteins:
[1] "P01624" "P0DOY3"
```

The exact numbers may differ depending on the protein file and protein list used.
