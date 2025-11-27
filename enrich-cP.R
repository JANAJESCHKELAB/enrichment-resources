#!/usr/bin/env Rscript
# enrich-cP.R
# By: Zian Liu
# Title: simple enrichment analyses via clusterProfiler
# Last updated: 2025-11-27

## TODO
## Add more visualization options that match real life use cases
## KEGG is not well-tested yet
## Add enrichr() universal function


# Add arguments
require(optparse)
opt.list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
    help="Input gene list file", metavar="character"
  ),
  make_option(c("-j", "--ontology"), type="character", default=NULL,
    help="Input ontology file (such as .gmt from Broad Institute)", metavar="character"
  ),
  make_option(c("-k", "--gosubcategory"), type="character", default="CC",
    help="For GO: sub-category (default=%default; can also be BP or MF)", metavar="character"
  ),
  make_option(c("-o", "--out"), type="character", default=NULL,
    help="Output directory/prefix", metavar="character"
  ),
  make_option(c("-e", "--enrichment"), type="character", default="GO",
    help="Type of enrichment test (default=%default; currently supports GO, GSEA)",
    metavar="character"
  ),
  make_option(c("-g", "--genenaming"), type="character", default='SYMBOL',
    help="The type of input gene naming, ALL CAP (default=%default)",
    metavar="character"
  ),
  make_option(c("-s", "--species"), type="character", default='human',
    help="Species (default=%default; only 'human' and 'mouse' supported for now)",
    metavar="character"
  ),
  make_option(c("-p", "--pvalue"), type="numeric", default=0.01,
    help="Maximum p-value cutoff (default=%default)"
  ),
  make_option(c("-q", "--qvalue"), type="numeric", default=0.05,
    help="Maximum q-value (adjusted p) cutoff (default=%default)"
  ),
  make_option(c("-V", "--visualize"), type="logical", default=TRUE,
    help="Whether to visualize results (default=%default)"
  )
)
opt.parser = OptionParser(option_list=opt.list);
opt = parse_args(opt.parser);

# Load packages & set up some environmental variables
require(clusterProfiler)
# Load species-specific databases
if (opt$species == 'human') {
  require(org.Hs.eg.db)
  data(geneList, package="DOSE")
  orgdb <- org.Hs.eg.db
} else if (opt$species == 'mouse') {
  require(org.Mm.eg.db)
  data(geneList, package="DOSE")
  orgdb <- org.Mm.eg.db
}
require(ggplot2)
require(ggarchery)
require(aplot)
require(enrichplot)
out.enrich <- paste0(opt$out, 'enrichment-results.csv')
if (opt$visualize) {
  out.plot <- opt$out
} else {
  out.plot <- NULL
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

readGenes <- function(filename, orgdb, enrichtype='GO', genenametype='SYMBOL') {
  # Load genes. Should be straightforward
  # If enrichtype is GSEA, load dataframe instead. Note that this has very
  # specific formatting requirement
  genes.df <- read.csv(filename, header=FALSE)
  if (enrichtype == 'GSEA') {
    # Col1 is some type of gene ID; Col2 is numeric
    genes.df[,1] <- as.character(genes.df[,1])
    names(genes.df) <- c(genenametype, "value")
    # Convert gene names
    if (genenametype != 'ENTREZID') {
      genes.converted <- convertGenenames(
        as.character(genes.df[,1]), genenametype, orgdb
      )
      genes.df <- dplyr::left_join(genes.df, genes.converted)
    }
    genes.list <- genes.df$value
    names(genes.list) = genes.df$ENTREZID
    genes.list = sort(genes.list, decreasing = TRUE)
  } else if (enrichtype %in% c('GO', 'KEGG')) {
    #genes.list <- genes.df[["V1"]]
    genes.list <- convertGenenames(genes.df[["V1"]], genenametype, orgdb)
  }
  return (genes.list)
}

runEnrichment <- function(
  geneobj, enrichtype='GO', globalgenes=NULL, orgdb=NULL, ontology='CC',
  file.ontology=NULL, species=NULL,
  pcutoff=0.01, qcutoff=0.05, filename=NULL) {
  # Enrichment wrapper
  # geneobj is output from previous (either dataframe or named list)
  # enrichment type should be c('GO', 'GSEA', 'KEGG', ...)
  # Some variables are type specific
  # `filename` if you want to save results
  if (enrichtype == 'GO') {
    res.enrichment <- enrichGO(
      gene = geneobj$ENTREZID,
      universe = names(globalgenes),
      OrgDb = orgdb,
      ont = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff = pcutoff,
      qvalueCutoff = qcutoff,
      readable = TRUE
    )
  } else if (enrichtype == 'KEGG') {
    if (species == 'human') {
      kegg.key = 'hsa'
    } else if (species == 'mouse') {
      kegg.key = 'mmu'
    }
    res.enrichment <- enrichKEGG(
      gene = geneobj$ENTREZID,
      organism = kegg.key,
      keyType = 'kegg',
      pvalueCutoff = pcutoff,
      universe = names(globalgenes)
    )
  } else if (enrichtype == 'GSEA') {
    gmt.df <- read.gmt(file.ontology)
    res.enrichment <- GSEA(
      geneobj,
      exponent = 1,
      minGSSize = 10,
      maxGSSize = 500,
      eps = 1e-10,
      pvalueCutoff = pcutoff,
      pAdjustMethod = "BH",
      TERM2GENE = gmt.df,
      TERM2NAME = NA,
      verbose = TRUE,
      seed = FALSE,
      by = "fgsea"
    )
  }
  if (!is.null(filename)) {
    write.csv(res.enrichment, filename)
  }
  return(res.enrichment)
}

visualizeResults <- function(res.enrich, enrichtype='GO', fileprefix=NULL) {
  # Basic visualization for enrichment results.
  # Write `filename` if you want to save it. Would be either name or prefix
  # depending on the plots generated.
  if (enrichtype %in% c('GO', 'KEGG')) {
    if (enrichtype == 'GO') {
      goplot(res.enrich)
      if (!is.null(fileprefix)) {
        ggsave(filename=paste0(fileprefix, enrichtype, 'plot.pdf'))
      }
    }
    barplot(res.enrich, showCategory=20)
    if (!is.null(fileprefix)) {
      ggsave(filename=paste0(fileprefix, enrichtype, 'Barplot.pdf'))
    }
  } else if (enrichtype == 'GSEA') {
    ridgeplot(res.enrich)
    p1 <- gseaplot(res.enrich, geneSetID = 1, by = "runningScore", title = res.enrich$Description[1])
    p2 <- gseaplot(res.enrich, geneSetID = 1, by = "preranked", title = res.enrich$Description[1])
    p3 <- gseaplot(res.enrich, geneSetID = 1, title = res.enrich$Description[1])
    plot_list(p1, p2, p3, ncol=1, tag_levels='A')
    if (!is.null(fileprefix)) {
      ggsave(filename=paste0(fileprefix, enrichtype, 'Ridgeplot.pdf'))
    }
  }
  # These are shared between GO/GSEA
  upsetplot(res.enrich)
  if (!is.null(fileprefix)) {
    ggsave(filename=paste0(fileprefix, enrichtype, 'Upsetplot.pdf'))
  }
}


## Main script below

gene.df <- readGenes(
  opt$input, orgdb, opt$enrichment, opt$genenaming
)
eres <- runEnrichment(
  geneobj=gene.df,
  enrichtype=opt$enrichment,
  globalgenes=geneList,
  file.ontology=opt$ontology,
  orgdb=orgdb,
  ontology=opt$gosubcategory,
  species=opt$species,
  pcutoff=opt$pvalue,
  qcutoff=opt$qvalue,
  filename=out.enrich
)
visualizeResults(
  eres, enrichtype=opt$enrichment, fileprefix=out.plot
)
