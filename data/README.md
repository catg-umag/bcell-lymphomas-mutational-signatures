# Data

This directory contains the minimal data of the project for running the notebooks.

## snv_list.csv
Each SNV used in the analyses. Reference: GRCh38.

### Fields
- `sample`: sample ID
- `lymph`: lymphoma type
- `sublymph`: lymphome subtype (mutated or unmutated for CLL)
- `chrom`: chromosome
- `pos`: position
- `ref`: reference nucleotide
- `alt`: alternative nucleotide
- `substitution`: subsitution in format `(ref)>(alt)`, converted to pyrimidines
- `context`: context (1 nucleotide each side)
- `context2`: context (2 nucleotides each side)
- `tumor_vaf`: VAF of tumor
- `gene`: gene(s) present in the mutation
- `ig`: mutation on an immunoglobulin gene
- `location`: region type
- `effect`: SNV effect


## signatures_nmf.Rdata
Pre-calculated signatures.