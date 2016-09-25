# Set of R scripts for analyzing math scores
# Scott Ness, Sep 2016 sness@salud.unm.edu
#

#### Start analysis loop

sample.name = samplenames[s]
sample.name

datadir = paste(outputdir, sample.name, "/", sep="")
datadir
mkdirs(datadir)

tumor_file = filenames[s]
tumor_file

## Perform Pre-Filtering of VCF Files

## Pre-filter tumor file
## Define prefilters
isSnp <- function(x) {
  refSnp <- nchar(ref(x)) == 1L
  a <- alt(x)
  altSnp <- elementLengths(a) == 1L
  ai <- unlist(a[altSnp]) # all length 1, so unlisting is 1:1 map
  altSnp[altSnp] <- nchar(ai) == 1L & (ai %in% c("A", "C", "G", "T"))
  refSnp & altSnp
}
minReadDepth <-function(x) {
  dp <- geno(x)$DP
  test <- (dp > tumor_min_read_depth)
  !is.na(test) & test
}
simpleGenotype <- function(x) {
  GT <- geno(x)$GT
  GTtest <- (GT %in% c("0/0", "0/1", "1/1"))
  !is.na(GTtest) & GTtest
}
altBothStrands <- function(x) {
  altfor <- geno(x)$SAF
  altrev <- geno(x)$SAR
  test <- (altfor > tumor_min_for_rev_reads) & (altrev > tumor_min_for_rev_reads)
  !is.na(test) & test
}
alleleFreq <- function(x) {
  af <- as.numeric(geno(x)$AF)
  if (tumor_vcf_old == TRUE) {
    af <- (as.numeric(geno(x)$AO) / as.numeric(geno(x)$RO))
  }
  test <- (af > tumor_min_var_allele_freq)
  !is.na(test) & test
}
rules <- list(isSnp=isSnp, simpleGenotype=simpleGenotype)
if (tumor_require_min_read_depth == TRUE) {
  rules <- list(isSnp=isSnp, simpleGenotype=simpleGenotype, minReadDepth=minReadDepth)
}
if (tumor_require_alt_reads_both_strands == TRUE) {
  rules <- list(isSnp=isSnp, simpleGenotype=simpleGenotype, altBothStrands=altBothStrands)
}
if (tumor_require_min_var_allele_freq == TRUE) {
  rules <- list(isSnp=isSnp, simpleGenotype=simpleGenotype, alleleFreq=alleleFreq)
}
if (tumor_require_min_read_depth == TRUE & tumor_require_alt_reads_both_strands == TRUE) {
  rules <- list(isSnp=isSnp, simpleGenotype=simpleGenotype, minReadDepth=minReadDepth, altBothStrands=altBothStrands)
}
if (tumor_require_min_read_depth == TRUE & tumor_require_min_var_allele_freq == TRUE) {
  rules <- list(isSnp=isSnp, simpleGenotype=simpleGenotype, minReadDepth=minReadDepth, alleleFreq=alleleFreq)
}
if (tumor_require_alt_reads_both_strands == TRUE & tumor_require_min_var_allele_freq == TRUE) {
  rules <- list(isSnp=isSnp, simpleGenotype=simpleGenotype, minReadDepth=minReadDepth, alleleFreq=alleleFreq)
}
if (tumor_require_min_read_depth == TRUE & tumor_require_alt_reads_both_strands == TRUE & tumor_require_min_var_allele_freq == TRUE) {
  rules <- list(isSnp=isSnp, simpleGenotype=simpleGenotype, minReadDepth=minReadDepth, altBothStrands=altBothStrands, alleleFreq=alleleFreq)
}

rules
filters <- FilterRules(rules)

## Apply filters
tumor_destination.file <- tempfile()
tabix.file <- TabixFile(tumor_file, yieldSize=10000)
filterVcf(tabix.file, "hg19", tumor_destination.file, filters=filters, verbose=TRUE)
#filterVcf(tabix.file, "b37", tumor_destination.file, filters=filters, verbose=TRUE)

## Read filtered tumor VCF file
tumor_all.vcf <- readVcf(tumor_destination.file, "hg19")
seqlevelsStyle(tumor_all.vcf) <- "UCSC"
seqlevels(tumor_all.vcf, force=TRUE) <- seqlevels(tumor_all.vcf)[seqlevels(tumor_all.vcf) != "chrM"]   # remove chrM

tumor_all.vcf
tumor.no.variants <- length(rownames(tumor_all.vcf))
tumor.no.variants
stopifnot(tumor.no.variants > 1)

# Collect header information from tumor VCF file
hdr <- scanVcfHeader(tumor_file)
file_format_tumor <- meta(hdr)$META["fileformat",]
file_source_tumor <- meta(hdr)$META["source",]

# Save tumor vcf file
outputfile = paste(datadir, sample.name, "_tumor_all.vcf", sep="")
writeVcf(tumor_all.vcf, outputfile, index = FALSE)


## Calculate Mutant-Allele Tumor Heterogeneity (MATH) score
if (use_all_variants_for_MATH == TRUE) {
  GT.tumor.unique <-geno(tumor_all.vcf)$GT  # Get genotype calls for all variants
  tumor.hetero <- GT.tumor.unique == "0/1"     # Limit the analysis only to the heterozygous calls
  tumor_all.vcf.hetero <- tumor_all.vcf[tumor.hetero]
} else {
  GT.tumor.unique <-geno(tumor_only.vcf)$GT  # Get genotype calls for tumor-specific variants
  tumor.hetero <- GT.tumor.unique == "0/1"     # Limit the analysis only to the heterozygous calls
  tumor_all.vcf.hetero <- tumor_only.vcf[tumor.hetero]
}
#math.vaf <- unlist(info(tumor_all.vcf.hetero)$AF)
if (vcfs_from_GATK == TRUE) {
  math.vaf <- sapply(geno(tumor_all.vcf.hetero)$AD, function(x) x[2]/(x[1]+x[2]))
} else {
  math.vaf <- round(as.numeric(geno(tumor_all.vcf.hetero)$AF), digits = 2)
}

head(math.vaf)
vaf2 <- (math.vaf > MATH_min_vaf) & (math.vaf < MATH_max_vaf)
vaf2 <- math.vaf[vaf2]
vaf2_length <- length(vaf2)
vaf_mad <- mad(vaf2)
vaf_med <- median(vaf2)
MATH_score <- (100 * vaf_mad)/vaf_med
MATH_score

sample.results <- data.frame(SAMPLES=sample.name, MATH=MATH_score, VARIANTS=tumor.no.variants)
sample.results

ms_data[s,] <- sample.results


#### ANALYZE TUMOR NONCODING VARIANTS
tumor_only.vcf <- tumor_all.vcf.hetero

## Annotate the tumor-specific variants
# Locate variants with respect to genes and coding regions
noncoding <- unique(locateVariants(tumor_only.vcf, txdb, AllVariants()))

noncoding.summary <- summary(noncoding$LOCATION)

no.splice <- as.matrix(noncoding.summary)[1]
no.intron <- as.matrix(noncoding.summary)[2]
no.five <- as.matrix(noncoding.summary)[3]
no.three <- as.matrix(noncoding.summary)[4]
no.coding <- as.matrix(noncoding.summary)[5]
no.intergenic <- as.matrix(noncoding.summary)[6]
no.promoter <- as.matrix(noncoding.summary)[7]

no.noncoding <- length(seqnames(noncoding))

no.splice
no.intron
no.five
no.three
no.coding
no.intergenic
no.promoter
no.noncoding

if (vcfs_from_GATK == TRUE) {
    vafs <- sapply(geno(tumor_only.vcf)$AD, function(x) x[2]/(x[1]+x[2]))
} else {
    vafs <- round(as.numeric(geno(tumor_only.vcf)$AF), digits = 2)
}
tumor_only.vcf.data <- data.frame(CDSLOC.names=rownames(tumor_only.vcf), SNP=paste(seqnames(tumor_only.vcf), start(tumor_only.vcf), sep="_"), VAF=round(vafs, digits = 2), DP=geno(tumor_only.vcf)$DP)
rownames(tumor_only.vcf.data) <- NULL
colnames(tumor_only.vcf.data) <- c("CDSLOC.names", "SNP", "VAF", "DP")
head(tumor_only.vcf.data)

# Add CDSLOC to files
snp_locs <- tumor_only.vcf.data$SNP
noncoding_locs <- paste(seqnames(noncoding), start(noncoding), sep="_")
m <- match(noncoding_locs, snp_locs)
snp_names <- tumor_only.vcf.data$CDSLOC.names[m]
mcols(noncoding)$CDSLOC.names <- snp_names

# Add VAF and DP
snp_vafs <- tumor_only.vcf.data$VAF[m]
snp_DPs <- tumor_only.vcf.data$DP[m]
mcols(noncoding)$VAF <- unlist(as.matrix(snp_vafs))
mcols(noncoding)$DP <- as.numeric(snp_DPs)

# Add gene symbols
m <- match(mcols(noncoding)$GENEID, egSYMBOL$gene_id)
symbols <- egSYMBOL$symbol[m]
mcols(noncoding)$SYMBOL <- symbols

# Remove troubling columns
mcols(noncoding)$PRECEDEID <- "<NA>"
mcols(noncoding)$FOLLOWID <- "<NA>"

# Add QUAL and FILTER fields from original vcf file
noncoding.ov <- findOverlaps(noncoding, rowRanges(tumor_only.vcf))
hits <- as.data.frame(noncoding.ov)
noncoding.ov.vcf <- tumor_only.vcf[hits[,2]]
quals <- rowRanges(noncoding.ov.vcf)$QUAL
filts <- rowRanges(noncoding.ov.vcf)$FILTER

mcols(noncoding)$QUAL <- quals
mcols(noncoding)$FILTER <- filts

## Summarize the number of variants by gene:
idx <- sapply(split(mcols(noncoding)$QUERYID, mcols(noncoding)$SYMBOL), unique)
genecounts <- sapply(idx, length)

outputfile = paste(datadir, "Noncoding_genecounts.txt", sep="")
write.table(genecounts, file=outputfile, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

noncoding.data.frame <- as.data.frame(noncoding, row.names = NULL, optional = FALSE)  # Use later for noncoding tumor mutations report


# All tumor gene mutations report
keep <- !is.na(noncoding.data.frame$GENEID)
#keep <- !noncoding.data.frame$SYMBOL == "<NA>"
xx <- summary(keep)
if ("TRUE" %in% names(xx)) {
    noncoding.data.frame.genes <- noncoding.data.frame[keep,]
    
    # Sort by SYMBOL
    x <- order(noncoding.data.frame.genes$SYMBOL, decreasing = FALSE)
    noncoding.data.frame.genes <- noncoding.data.frame.genes[x,]
    
    outputfile <- paste(datadir, "tumor_all_mutations.txt", sep="")
    write.table(noncoding.data.frame.genes, file=outputfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    
    tumor_all_coding_report = TRUE
} else {
    tumor_all_coding_report = FALSE
}


# End of analysis

