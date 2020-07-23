# BENEATH THE REPO (kind of READ ME)

### Prepearing working environment:
source("csourse.R")
source("cutils.R")
work.dir <- "~/Desktop/"
dir.create(work.dir, showWarnings = F)


### Data preprocessing:
# *
# Load data obtained like in ...R (Formal requirments to data structure and properties: ...)
load("~/Documents/immGen/final_objects/243_es.top12k.Rda")
# *
Biobase::pData(es.top12k)$sample <- colnames(es.top12k)
# If you don't have entrez IDs in your annotation, take gene symbols
entrez <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
                                column = "ENTREZID",
                                keytype = "SYMBOL",
                                keys=as.character(Biobase::fData(es.top12k)$gene))
# *
# Extracting metabolic genes
exprsMetab <- extractMetabolicGenes(exprs = Biobase::exprs(es.top12k),
                                    entrezIDs = entrez)
# *
#' Get proposed parameters for the input data.
#' It can take a long time. If you want to skip this step you can use k = 64 and base = 0.4.
#'
#' @param exprs Matrix with expression.
#' @return k is the number
#' @return base is the number
#' @example findParameters(exprs = Biobase::exprs(es.top12k))
proposedParameters <- findParameters(exprsMetab = exprsMetab)
proposedParameters$k
proposedParameters$base



### Initial clustering:

#' Initial clustering. Sets the number of initital clusters.
#'
#' @param es value: ExpressionSet
#' @param organism value: choose either `mouse`, either `human`
#' @param repeats e.g. you may use gsub() ti collapce replicas or use any custon annotation
#' @param initialNumber
#' @param showInitialClustering show heatmap of initial clusters
#' @return The table sum of \code{x} and \code{y}.
#'
#' # Currently availble for Mus musculus and GAM method only
preCluster <- preClustering(es = es.top12k,
                            repeats = Biobase::pData(es.top12k)$sample,
                            entrezIDs = Biobase::fData(es.top12k)$entrez,
                            initialNumber = 5, # proposedParameters$k
                            showInitialClustering = T)
str(preCluster)
str(preCluster$network)
str(preCluster$gene.exprs_orig)
str(preCluster$gene.exprs)
str(preCluster$curCenters)

library(GAM) # get.vertex.attributes
library(plyr) # rename (nb! not dplyr)

#' GAM-clustering per se
#'
#' @param gene.exprs value: ExpressionSet
#' @param base reflects the degree of retention of edges in the active module
#' @return The table sum of \code{x} and \code{y}.
gamCluster <- gamClustering(network = preCluster$network,
                            gene.exprs = preCluster$gene.exprs,
                            curCenters = preCluster$curCenters,
                            repeats = Biobase::pData(es.top12k)$sample,
                            base = 0.4, # proposedParameters$base
                            work.dir = work.dir,
                            batch.script = "/home/octopus/R-studio/nclust/sgmwcs-slurm-batch-0.9.5",
                            showIntermediateClustering = TRUE,
                            saveSession = TRUE)
gamCluster$k
gamCluster$revs
gamCluster$curRev
gamCluster$curRev$modules
session::restore.session(file=paste0(work.dir, "/session.RDa"))


### Visualizing gam-clustering results:

library(parallel)
#' Get files with modules' graphs.
#'
#' @param curRev value:
#' @return results of these functions calling can be seen in work.dir
getGraphs(curRev = gamCluster$curRev,
          work.dir = work.dir)

#' Get files with modules' heatmaps.
#'
#' @param curRev value:
#' @return results of these functions calling can be seen in work.dir
getHeatmaps(curRev = gamCluster$curRev,
            gene.exprs = preCluster$gene.exprs,
            work.dir = work.dir)

#' Get files with modules' gene tables.
#'
#' @param curRev value:
#' @return results of these functions calling can be seen in work.dir
getGeneTables(curRev = gamCluster$curRev,
              nets = gamCluster$nets,
              gene.exprs = preCluster$gene.exprs,
              work.dir = work.dir,
              organism = "mouse")

#' Get files with modules' annotations by canonical pathways.
#'
#' @param modules The modules derived from @gamClustering function calling.
#' @param files The files derives from @getGeneTables function calling.
#' @return The files with pathways...
#' @example annotateModules(modules, files)
annotateModules(modules)



