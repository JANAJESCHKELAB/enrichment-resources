#!/usr/bin/env Rscript
# enrich-cP.R
# By: Zian Liu
# Title: Basic enrichment analyses via clusterProfiler
# Last updated: 2025-11-27


# Add arguments
require(optparse)
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
    help="Input gene list file", metavar="character"
  ),
  make_option(c("-o", "--out"), type="character", default=NULL,
    help="Output directory/prefix", metavar="character"
  ),
  make_option(c("-e", "--enrichment"), type="character", default="GO",
    help="Type of enrichment test (default=%default; only 'GO' supported for now)",
    metavar="character"
  ),
  make_option(c("-s", "--species"), type="character", default='human',
    help="Species (default=%default; only 'human' supported for now)](",
    metavar="character"
  ),
  make_option(c("-p", "--pvalue"), type="numeric", default=0.01,
    help="Maximum p-value cutoff (default=%default)"
  ),
  make_option(c("-q", "--qvalue"), type="numeric", default=0.05,
    help="Maximum q-value cutoff (default=%default)"
  ),
  make_option(c("-V", "--visualize"), type="logical", default=TRUE,
    help="Whether to visualize results (default=%default)"
  )
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Make input variables
input.genes <- opt$input
type.enrichment <- opt$enrichment  # Add GSEA later
species <- opt$species  # TODO: add more options?
out.prefix <- opt$out
pcutoff <- opt$pvalue
qcutoff <- opt$qvalue

require(clusterProfiler)
# Load species-specific databases
if (species == 'human') {
  require(org.Hs.eg.db)
  data(geneList, package="DOSE")
  orgdb <- org.Hs.eg.db
}
require(ggplot2)
require(ggarchery)

# These are generated automatically
out.enrich <- paste0(out.prefix, 'enrichment-results.csv')
if (opt.visualize) {
  out.plot <- paste0(out.prefix, 'enrichment.pdf')
} else {
  out.plot <- NULL
}


readGenes <- function(filename) {
  # Load genes. Should be straightforward
  genes.df <- read.delim(filename, header=FALSE)
  return (genes.df[["V1"]])
}

## TODO: not done yet
readGenesGSEA <- function(filename) {
  # For GSEA. Note that this has very specific formatting requirements
  genes.df <- read.csv(filename, header=FALSE)
  ## TODO

  ## /TODO
}

convertGenenames <- function(genevector, keytype, orgdb) {
  # Convert to other formats.
  # input_type has to be a supported format
  typeall <- c(
    'ACCNUM', 'ALIAS', 'ENSEMBL', 'ENSEMBLPROT', 'ENSEMBLTRANS', 'ENTREZID', 'ENZYME', 'EVIDENCE',
    'EVIDENCEALL', 'GENENAME', 'GENETYPE', 'GO', 'GOALL', 'IPI', 'MAP', 'OMIM', 'ONTOLOGY',
    'ONTOLOGYALL', 'PATH', 'PFAM', 'PMID', 'PROSITE', 'REFSEQ', 'SYMBOL', 'UCSCKG', 'UNIPROT'
  )
  if (!keytype %in% typeall) {
    stop("The input type is not an clusterProfiler accepted type!")
  }
  translate.df <- bitr(
    genevector, fromType = keytype,
    toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
    OrgDb = orgdb
  )
  return (translate.df)
}

runEnrichment <- function(
  gene, enrichtype='GO', globalgenes=NULL, orgdb=NULL, ontology='CC',
  pcutoff=0.01, qcutoff=0.05, filename=NULL) {
  # Enrichment wrapper
  # enrichment type should be c('GO', 'GSEA', ...)
  # Some variables are type specific
  # `filename` if you want to save results
  ## TODO: This is WIP

  if (enrichtype == 'GO') {
    res_enrichment <- enrichGO(
      gene = gene,
      universe = names(globalgenes),
      OrgDb = orgdb,
      ont = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff = pcutoff,
      qvalueCutoff = qcutoff,
      readable = TRUE
    )
  }
  if (!is.null(filename)) {
    write.csv(res_enrichment, filename)
  }
  return(res_enrichment)
}

visualizeResults <- function(res_enrich, filename=NULL) {
  # Basic visualization for enrichment results.
  # Write `filename` if you want to save it. Would be either name or prefix
  # depending on the plots generated.
  goplot(res_enrich)
  if (!is.null(filename)) {
    ggsave(filename=filename)
  }
}

genes <- readGenes(input.genes)
gene.df <- convertGenenames(genes, "SYMBOL", orgdb)
eres <- runEnrichment(
  gene=gene.df$ENTREZID,
  enrichtype=type.enrichment,
  globalgenes=geneList,
  orgdb=orgdb,
  pcutoff=pcutoff,
  qcutoff=qcutoff,
  filename=out.enrich
)
visualizeResults(eres, out.plot)
