# Introduction

This is an effort to gather some easy-to-use gene
enrichment resources for the entire ACEE.
`example-*/` contain exemplar gene sets and anticipated
results to test the script with.

If you find this useful, please help contribute.


## Using the script

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


## External resources

* [Gene Ontology](https://geneontology.org/):
They have a built-in GO analysis software and some visualization
tools.
* [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp): 
if you want to run GSEA you should download one of these gene sets
first. Need to enter personal information.
* [PAN-GO](https://functionome.geneontology.org/):
"all annotated functional characteristics for human protein-coding genes".
They have some online tools too.



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
