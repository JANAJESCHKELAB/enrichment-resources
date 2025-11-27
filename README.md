# Introduction

This is an effort to make a somewhat streamlined and
easy-to-use enrichment script for ACEE.
`example-*/` contain exemplar gene sets and anticipated
results to test the script with.


## Use

`./enrich-cP.R <options>`

Use `--help` to check available options.

If you just want to run GO,
`./enrich-cP.R <inputfile> <output_directory>` is sufficient.


## Current features/points of notice

* Supports good-old-fashioned over-representation analysis (ORA)
with Gene Ontology (GO) terms
* Supports GSEA analysis; make sure to download a .gmt gene set file
from [GSEA](https://www.gsea-msigdb.org) as they are needed
* Only supports human or mouse for now


## TODO

* Add support for custom databases (such as using local PAN-GO)


## References

* [Gene Ontology](https://geneontology.org/)
* [clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler)
* [GSEA](https://www.gsea-msigdb.org)
