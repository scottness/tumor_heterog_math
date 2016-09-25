# math_score_functions.R
#
# Scott Ness, Sep 2016 sness@salud.unm.edu

## Test for output directory
# from CRAN
library(R.utils)
library(checkmate)

# Set up the output directory and file names
analysis.date <- format(Sys.Date(), "%y%m%d")
workingdir <- getwd()

#### Load packages and functions
# from CRAN
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(plot3D)
library(rgl)
library(xtable)
library(diagram)
library(lattice)

# from Bioconductor
library(edgeR)
library(DESeq)
library(genefilter)
library(VariantAnnotation)
library(BSgenome)
library(AnnotationHub)
library(pathview)
library(gage)
library(topGO)
library(biomaRt)
library(org.Hs.eg.db)
egSYMBOL <- toTable(org.Hs.egSYMBOL)
egREFSEQ <- toTable(org.Hs.egREFSEQ)
egENSEMBL <- toTable(org.Hs.egENSEMBL)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
if (!exists("txdb")) {
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
}

# 
