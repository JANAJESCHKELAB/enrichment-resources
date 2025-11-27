# Introduction

This is an effort to make a somewhat streamlined and
easy-to-use enrichment script for ACEE.
`example-*/` contain exemplar gene sets and anticipated
results to test the script with.


## Use

First, prepare your inputs:
* Input gene list should be a newline-separated
list of genes in a plain text file. Any clusterProfiler-supported
gene ids (symbols, Entrez IDs, ...) are fine; default is gene symbols.
* Specifically for GSEA, input should be an unnamed (no column/row names) 
.csv file, with column 1 being the gene ids (same requirements as above),
and column 2 being "values" (whatever that may be in your case).
* If you are using an enrichment method that requires an external
downloaded ontology, make sure to download it. This includes GSEA
(requires .gmt files from Broad Institute).

Then, run:

`./enrich-cP.R <options>`

Use `--help` to check available options.

If you just want to run GO,
`./enrich-cP.R <input_genelist> <output_directory>` is sufficient.


## Current features/points of notice

* Supports good-old-fashioned over-representation analysis (ORA)
with Gene Ontology (GO) terms
* Supports GSEA analysis; make sure to download a .gmt gene set file
from [GSEA](https://www.gsea-msigdb.org)
* Only supports human or mouse for now


## TODO

* Add support for custom databases (such as using local PAN-GO)


## References

* [Gene Ontology](https://geneontology.org/)
* [clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler)
* [GSEA](https://www.gsea-msigdb.org)
