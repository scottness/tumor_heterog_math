# Set of R scripts for analyzing math scores
# Scott Ness, Sep 2016 sness@salud.unm.edu
# 
#### Description
# 

#### Update the information in the next few lines, then run...
#
## Define the file locations and analysis parameters
analysis.name = "math_analysis"      # Define this for each project or set of samples
path_to_save = "PROVIDE A PATH TO SAVE THE OUTPUT"   

# Location, files to be analyzed and parameters
path_to_files = "PROVIDE A PATH TO THE DATA FILES"
vcfs_from_GATK = FALSE               # TRUE if analyzing RNA-seq VCF files
recursive_value = FALSE              # TRUE to look through sub-folders
pattern_to_match = ".vcf.gz$"   # Edit this to find all the appropriate VCF files

samples_file = "PROVIDE A PATH TO THE SAMPLE DESCRIPTION FILE"
# Samples data should be a CSV file with these fields:
# UID               (unique identifier, used as sample name)
# INCLUDE           (any samples not marked "yes" will be excluded)
# GROUP_LETTER      (Letter code, e.g.: A,B,C,D,E)
# GROUP_NO          (Number, e.g.: 1,2,3)
# PRIMARY_FILE      (file base name)
# Additional columns are allowed

include_marked_only = TRUE        # If true, only use samples with "yes" in the samples data INCLUDE column



####################
## These are advanced settings that rarely change
####################

#### Sequence analysis settings ####
## Define tumor filter parameters
filter_tumor_file = TRUE
tumor_vcf_old = FALSE
# Should be SNPs (no InDels or other)
# Should have simple genotype
tumor_require_min_read_depth = TRUE
tumor_min_read_depth = 50                     # Minimum read depth
tumor_require_alt_reads_both_strands = TRUE   # Should have alternate allele reads on both strands
tumor_min_for_rev_reads = 4                   # must be greater than this, setting to 3 means minimum 4 reads each direction
tumor_require_min_var_allele_freq = TRUE
tumor_min_var_allele_freq = 0.05              # use if only looking for germline variants
remove_low_dp_from_plots = FALSE              # cut off low dp variants for plots
use_vaf_to_define_germline = FALSE            # define germline hetero and homozygous by VAF (not for RNA-seq)

# MATH score calulation settings
use_all_variants_for_MATH = TRUE              # Default = TRUE.  FALSE excludes the germline variants from MATH score
MATH_min_vaf = 0.05                           # lower limit for VAF included in MATH score
MATH_max_vaf = 0.75                           # upper limit for VAF included in MATH score

## Values for PCA plot section
use_ranked_snps = FALSE            # Limit heatmap to only most different SNPs
min_samples_with_snp = 1          # At least this many samples must have the variant for it to appear in the heatmap
heatmap_min_vaf = 0.05            # Lower limit for variants to include in heatmap
heatmap_max_vaf = 0.9             # Upper limit for variants to include in heatmap


####################
## Nothing below this should need editing
####################

# Define output directory for filtered VCF files

## Extract list of files to analyze
filenames <- list.files(path=path_to_files, pattern = "*.vcf.gz$", recursive = recursive_value, ignore.case = FALSE, full.names=TRUE)
filenames
found <- grepl(pattern_to_match, filenames, ignore.case = TRUE)  # find e.g. normal or tumor files only
filenames <- filenames[found]
filenames
filenames_all <- filenames
basenames <- basename(filenames)
basenames

# Read a file of tumor stage or other meta data
samples.data <- read.csv(samples_file, header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
samples.data

# Keep only files marked include = yes
if (include_marked_only == TRUE) {
  inc <- samples.data$INCLUDE == "yes"
  samples.data <- samples.data[inc,]
}
samples.data

keep <- samples.data$PRIMARY_FILE %in% basenames
summary(keep)
samples.data <- samples.data[keep,]

keep <- basenames %in% samples.data$PRIMARY_FILE
summary(keep)

basenames <- basenames[keep]
filenames <- filenames[keep]

ff <- t(data.frame(sapply(basenames, strsplit, split='_')))
ff
samplenames <- as.matrix(paste(ff[,1], ff[,2], sep="_"))
samplenames

# set up data file
ms_data <- data.frame(SAMPLES=samplenames, MATH=0, VARIANTS=0)

basenames
filenames

no.samples <- length(samplenames)
no.samples

# make sure the order of files and meta data is correct
m <- match(basenames, samples.data$PRIMARY_FILE)
m
samples.data <- samples.data[m,]

labels <- samplenames
labels

# Load all the libraries and functions
source("R/math_score_functions.R", echo=TRUE)

outputdir = paste(path_to_save, "/", analysis.name, "_", analysis.date, "/", sep="")
outputdir

mkdirs(outputdir)                           # Creates output directory if it doesn't exist, old data gets over-written if it does
testDirectory(outputdir, "rw")              # Test if output directory exists and is writeable
stopifnot(testDirectory(outputdir, "rw"))   # stop script if output directory not writeable
outputdir


## Analysis loop
for (s in 1:no.samples) {
source("R/math_score_script.R", echo=TRUE)
}


# Save ms data file
outputfile = paste(outputdir, "ms_data.csv", sep="")
write.table(ms_data, file=outputfile, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE, qmethod = "double")

# Add ms data to samples.data file
samples.data$MATH <- ms_data$MATH
samples.data$VARIANTS <- ms_data$VARIANTS

# Save complete sample data file
outputfile = paste(outputdir, "step1_sample_data.csv", sep="")
write.table(samples.data, file=outputfile, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE, qmethod = "double")



####################
## Prepare MATH score boxplot
####################

samples.data

# Perform t.test of different groups
group2 <- samples.data$GROUP_NO == 2
group2 <- samples.data[group2,]
group2
no.group2 <- length(group2[,1])

group3 <- samples.data$GROUP_NO == 3
group3 <- samples.data[group3,]
group3
no.group3 <- length(group3[,1])


ms.stage.2v3.t.test <- t.test(as.numeric(as.matrix(group2$MATH)), as.numeric(as.matrix(group3$MATH)), var.equal = TRUE)
ms.stage.2v3.t.test


plot(samples.data$GROUP_NO, (samples.data$MATH))


# Make a boxplot
y = as.numeric(as.matrix(samples.data$MATH))
grp = samples.data$GROUP_NO
title = "Mutant-Allele Tumor cosmic (MATH) Score by Stage"
xlabel = "Tumor Stage"
ylabel = "MATH Score (100*MAD/median VAF at heterozygous tumor loci)"
xlab1 <- paste("Stage II (n=", no.group2, ")", sep="")
xlab2 <- paste("Stage III (n=", no.group3, ")", sep="")
grp_names = c(xlab1, xlab2)
p_val1 <- paste("pval=", round(ms.stage.2v3.t.test$p.val, digits = 3), sep="")
y_min <- round(min(samples.data$MATH)-1,0)
y_max <- round(max(samples.data$MATH)+4, 0)
ylimits <- c(y_min,y_max)

name = paste(outputdir, "MATH_score_boxplot.jpg",  sep = "", collapse = NULL)
jpeg(name, height = 700, width = 700)
boxplot(y ~ grp, main = title, xlab = xlabel, ylab = ylabel, ylim = ylimits, varwidth = TRUE, names = grp_names)
par(xpd=TRUE)
yrange<-par("usr")[3:4]
ypos1<-(yrange[2]+diff(yrange)/40) - 4
segments(1,ypos1,2,ypos1)
text(1.5,(1+ypos1+diff(yrange)/40),p_val1)
if (ms.stage.2v3.t.test$p.val < 0.05) {
    text(1.5,(ypos1+diff(yrange)/40),"*",cex=2)
}
par(xpd=FALSE)
dev.off()

name = paste(outputdir, "MATH_score_boxplot.pdf",  sep = "", collapse = NULL)
pdf(name, height = 7.5, width = 7.5)
boxplot(y ~ grp, main = title, xlab = xlabel, ylab = ylabel, ylim = ylimits, varwidth = TRUE, names = grp_names)
par(xpd=TRUE)
yrange<-par("usr")[3:4]
ypos1<-(yrange[2]+diff(yrange)/40) - 4
segments(1,ypos1,2,ypos1)
text(1.5,(1+ypos1+diff(yrange)/40),p_val1)
if (ms.stage.2v3.t.test$p.val < 0.05) {
    text(1.5,(ypos1+diff(yrange)/40),"*",cex=2)
}
par(xpd=FALSE)
dev.off()



####################
## Prepare SNP PCA plots
####################

filenames_all_tumor <- list.files(path=outputdir, pattern = "tumor_all_mutations.txt$", recursive = TRUE, ignore.case = FALSE, full.names=TRUE)

filenames_all_tumor
samples.names <- samplenames


# Read the results from each of the Nof1 files
All <- lapply(filenames_all_tumor,function(i){
    read.delim(i, row.names = NULL)
})

length(All[[1]]$seqnames)


## Remove entries with no symbol and duplicates
for (i in 1:no.samples) {
    if (length(All[[4]]$seqnames) > 0) {
        keep <- !is.na(All[[i]]$SYMBOL)
        All[[i]] <- All[[i]][keep,]
        x <- All[[i]]$CDSLOC.names
        All[[i]] <- All[[i]][!duplicated(x),]
        All[[i]] <- All[[i]][,1:length(colnames(All[[i]]))]
    }
}

length(All[[1]]$seqnames)

# Generate table of number of mutations in each sample
lengths_all = matrix(nrow=no.samples, ncol=1)
for (i in 1:no.samples) {
    lengths_all[i] <- length(All[[i]]$seqnames)
}
rownames(lengths_all) <- samplenames
colnames(lengths_all) <- c("NUMBER")
lengths_all

# Determine the set of overlapping SNP calls
found <- as.matrix(All[[1]]$CDSLOC.names)
length(found)

for (i in 1:no.samples) {
    test <- as.matrix(All[[i]]$CDSLOC.names)
    found <- intersect(found, test)
}
no.common.genotypes <- length(found)
no.common.genotypes

# Build master table of variants
master <- rbind(All[[1]])
for (i in 2: no.samples) {
    keep <- !All[[i]]$CDSLOC.names %in% master$CDSLOC.names
    uni <- All[[i]][keep,]
    master <- rbind(master, uni)
}

o <- order(master$seqnames, master$start, decreasing = FALSE)
master <- master[o,]
rownames(master) <- NULL
head(master)

P <- lapply(All, function(x){
    x <- paste(master$seqnames, master$start, sep="_") %in% paste(x$seqnames, x$start, sep="_")
})
master2 <- data.frame(master)
for (i in 1:no.samples) {
    master2 <- data.frame(master2, P[[i]])
}
cols <- colnames(master)
colnames(master2) <- c(cols, samplenames)

w = length(colnames(master2))-no.samples
w

test <- (master2[,(w+1):(w+no.samples)])
x <- apply(test, 2, function(x) as.numeric(x))
count <- rowSums(test)

master2 <- data.frame(master2, SNP_COUNT=count)

o <- order(master2$SYMBOL)
master2 <- master2[o,]

genecounts <- as.matrix(tapply(master2$SNP_COUNT, master2$SYMBOL, sum, na.rm = TRUE))
keep <- !is.na(genecounts)
genecounts <- genecounts[keep,]
genecounts
o <- order(names(genecounts))
genecounts <- genecounts[o]
genecounts

outputfile = paste(outputdir, "All_tumor_genecounts_summary.txt", sep="")
write.table(genecounts, file=outputfile, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

# Make barchart
x <- as.matrix(genecounts)
colnames(x) <- c("counts")
o <- order(x, decreasing = FALSE)
x <- x[o,]
l <- length(x)
maintitle = "All Tumor Genecounts"

## Save barcharts
name = paste(outputdir, "All_tumor_genecounts_all_barchart.jpg",  sep = "", collapse = NULL)
jpeg(name, height = 700, width = 700)
if (l>24) {
    barchart(x[(l-24):l], horizontal=TRUE, xlab = "Frequency", main=maintitle)
} else {
    barchart(x, horizontal=TRUE, xlab = "Frequency", main=maintitle)
}
dev.off()

name = paste(outputdir, "All_tumor_genecounts_all_barchart.pdf",  sep = "", collapse = NULL)
pdf(name, height = 7.5, width = 7.5)
if (l>24) {
    barchart(x[(l-24):l], horizontal=TRUE, xlab = "Frequency", main=maintitle)
} else {
    barchart(x, horizontal=TRUE, xlab = "Frequency", main=maintitle)
}
dev.off()

# Save all tumor mutations data
outputfile = paste(outputdir, "All_tumor_mutations_data.txt", sep="")
write.table(master2, file=outputfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)



# ## Select for response-specific SNPs
group1 <- samples.data$GROUP_Letter == "B"
group2 <- samples.data$GROUP_Letter == "A"
w <- length(colnames(master2))
w
no.samples
master2.data1 <- t(master2[,(w-no.samples):(w-1)])
master2.data1 <- master2.data1[group1,]
master2.data1 <- t(master2.data1)
head(master2.data1)
l1 <- length(colnames(master2.data1))
master2.data2 <- t(master2[,(w-no.samples):(w-1)])
master2.data2 <- master2.data2[group2,]
master2.data2 <- t(master2.data2)
head(master2.data2)
l2 <- length(colnames(master2.data2))

w = (length(colnames(master2))-no.samples-1)

test1 <- master2.data1
x <- apply(test1, 2, function(x) as.numeric(x))
count1 <- rowSums(test1)
master2_1 <- data.frame(master2[,1:w], master2.data1, SNP_COUNT=count1)
head(master2_1)

test2 <- master2.data2
x <- apply(test2, 2, function(x) as.numeric(x))
count2 <- rowSums(test2)
master2_2 <- data.frame(master2[,1:w], master2.data2, SNP_COUNT=count2)
head(master2_2)

ranks <- data.frame(COUNT1=master2_1$SNP_COUNT, COUNT2=master2_2$SNP_COUNT)
ranks <- data.frame(ranks, RATE1=ranks$COUNT1/l1, RATE2=ranks$COUNT2/l2, DIFF=(ranks$COUNT1/l1)-(ranks$COUNT2/l2))
rownames(ranks) <- rownames(master2_1)
head(ranks)

hist(ranks$DIFF)

diff.low <- ranks$DIFF < -0.25
summary(diff.low)

diff.high <- ranks$DIFF > 0.35
summary(diff.high)

master2_low <- master2[diff.low,]
master2_high <- master2[diff.high,]

master2_ends <-rbind(master2_low, master2_high)
head(master2_ends)

o <- order(as.character(master2_ends$SYMBOL))
master2_ends <- master2_ends[o,]

master2b <- master2

if (use_ranked_snps == TRUE) {
    master2b <- master2_ends
}



# Build unified table of vafs
head(master2b)
# rownames(master)
keep <- master2b$SNP_COUNT >= min_samples_with_snp
summary(keep)
snp_list <- master2b[keep,]
snp_list <- data.frame(snp_list, INX=paste(snp_list$seqnames, "_", snp_list$start, sep=""))

uni_snp_list <- snp_list[!duplicated(snp_list$INX),]
o <- order(uni_snp_list$seqnames, uni_snp_list$start, decreasing = FALSE)
uni_snp_list <- uni_snp_list[o,]

# Build VAF lists
vafs <- NULL
for (i in 1:no.samples) {
    vafs[[i]] <- data.frame(INX=paste(All[[i]]$seqnames, "_", All[[i]]$start, sep=""), All[[i]]$VAF)
    colnames(vafs[[i]]) <- c("INX", samplenames[i])
    vafs[[i]] <- vafs[[i]][!duplicated(vafs[[i]]$INX),]
    keep <- vafs[[i]][,2] > 0.01
    vafs[[i]] <- vafs[[i]][keep,]
    keep <- vafs[[i]][,2] > heatmap_min_vaf
    vafs[[i]] <- vafs[[i]][keep,]
    keep <- vafs[[i]][,2] < heatmap_max_vaf
    vafs[[i]] <- vafs[[i]][keep,]
    both.names <- intersect(vafs[[i]]$INX, uni_snp_list$INX)
    diff.names <- unique(setdiff(uni_snp_list$INX, vafs[[i]]$INX))
    keep <- vafs[[i]]$INX %in% both.names
    summary(keep)
    vafs[[i]] <- vafs[[i]][keep,]
    colnames(vafs[[i]]) <- c("INX", "x")
    add_to1 <- data.frame(INX=diff.names, x=0)
    vafs[[i]] <- rbind(vafs[[i]], add_to1)
    rownames(vafs[[i]]) <- vafs[[i]]$INX
    o <- order(rownames(vafs[[i]]), decreasing=FALSE)
    vafs[[i]] <- vafs[[i]][o,]
    colnames(vafs[[i]]) <- c("INX", samplenames[i])
}

# Add VAF lists to VAF table
o <- order(uni_snp_list$INX, decreasing = FALSE)
uni_snp_list <- uni_snp_list[o,]
vaf_table <- data.frame(uni_snp_list[,1:15], SYMBOL=uni_snp_list$SYMBOL, INX=uni_snp_list$INX)
cols <- colnames(vaf_table)
for (i in 1:no.samples) {
    vaf_table <- data.frame(vaf_table, vafs[[i]][,2])
}
keep <- apply(vaf_table[,18:w],1,sum) > 0    # remove any rows with only zeros
vaf_table <- vaf_table[keep,]
colnames(vaf_table) <- c(cols, samplenames)
rownames(vaf_table) <- vaf_table$INX
head(vaf_table)



## Make PCA plots and build custom heatmap
w = length(colnames(vaf_table))
w
pdata2 <- as.matrix(vaf_table[,18:w])
rownames(pdata2) <- paste(vaf_table$CDSLOC.names, vaf_table$SYMBOL, sep="_")
head(pdata2)
IDs <- colnames(pdata2)
IDs
m <- match(IDs, samples.data$UID)
m
samples.data2 <- samples.data[m,]
samples.data

color.map <- function(group.col) { if (group.col=="1") "#dfdfdf" 
    else if (group.col=="2") "#000000" 
    else if (group.col=="3") "#FF0000"
    else if (group.col=="4") "#0000FF"
    else if (group.col=="5") "#FFCC00" 
    else "#FFFFFF" }
groupcolors <- unlist(lapply(samples.data$GROUP_NO, color.map))
groupcolors


# ## Perform PCA analysis (MDS)
name = paste(outputdir, "All_tumor_genecounts_VAF_MDS_plot.jpg",  sep = "", collapse = NULL)
jpeg(name, height = 700, width = 700)
plotMDS(pdata2, col = groupcolors)
dev.off()

name = paste(outputdir, "All_tumor_genecounts_VAF_MDS_dim12.pdf",  sep = "", collapse = NULL)
pdf(name, height = 7.5, width = 7.5)
mds1 <- plotMDS(pdata2, col = groupcolors)
dev.off()

name = paste(outputdir, "All_tumor_genecounts_VAF_MDS_dim23.pdf",  sep = "", collapse = NULL)
pdf(name, height = 7.5, width = 7.5)
mds2 <- plotMDS(pdata2, dim.plot=c(2,3), col = groupcolors)
dev.off()

name = paste(outputdir, "All_tumor_genecounts_VAF_MDS_dim13.pdf",  sep = "", collapse = NULL)
pdf(name, height = 7.5, width = 7.5)
mds3 <- plotMDS(pdata2, dim.plot=c(1,3), col = groupcolors)
dev.off()


# Make 3D MDS Plot
xdim <- mds1$x
ydim <- mds1$y
zdim <- mds2$y

scatter3D(xdim,ydim,zdim, cex = 2, col=c(groupcolors))

plot3d(xdim,ydim,zdim, cex = 2, size=20, col=c(groupcolors), bg="white")



# Save complete session info
sessionfile = paste(outputdir, "Session_Info.txt", sep="")
si <- capture.output(sessionInfo())
write.table(si, file = sessionfile, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)







####################
## Compare Additional Stage II samples
####################
rm(list=ls())


# Update the information in the next few lines, then run...
#
## Define the file locations and analysis parameters
analysis.name = "math_analysis_stageII"      # Define this for each project or set of samples
path_to_save = "PROVIDE A PATH TO SAVE THE RESULTS"   

# Location, files to be analyzed and parameters
path_to_files = "PROVIDE A PATH TO THE VCF FILES"
vcfs_from_GATK = FALSE          # TRUE if analyzing RNA-seq VCF files
recursive_value = FALSE         # TRUE to look through sub-folders
pattern_to_match = ".vcf.gz$"   # Edit to match all the VCF files

samples_file = "PROVIDE A PATH TO THE SAMPLE DESCRIPTION FILE"
# Samples data should be a CSV file with these fields:
# UID               (unique identifier, used as sample name)
# INCLUDE           (any samples not marked "yes" will be excluded)
# GROUP_LETTER      (Letter code, e.g.: A,B,C,D,E)
# GROUP_NO          (Number, e.g.: 1,2,3)
# PRIMARY_FILE      (file base name)
# Additional columns are allowed

include_marked_only = TRUE        # If true, only use samples with "yes" in the samples data INCLUDE column


####################
## These are advanced settings that rarely change
####################

#### Sequence analysis settings ####
## Define tumor filter parameters
filter_tumor_file = TRUE
tumor_vcf_old = FALSE
# Should be SNPs (no InDels or other)
# Should have simple genotype
tumor_require_min_read_depth = TRUE
tumor_min_read_depth = 50                     # Minimum read depth
tumor_require_alt_reads_both_strands = TRUE   # Should have alternate allele reads on both strands
tumor_min_for_rev_reads = 4                   # must be greater than this, setting to 3 means minimum 4 reads each direction
tumor_require_min_var_allele_freq = TRUE
tumor_min_var_allele_freq = 0.05              # use if only looking for germline variants
remove_low_dp_from_plots = FALSE              # cut off low dp variants for plots
use_vaf_to_define_germline = FALSE            # define germline hetero and homozygous by VAF (not for RNA-seq)

# MATH score calulation settings
use_all_variants_for_MATH = TRUE              # Default = TRUE.  FALSE excludes the germline variants from MATH score
MATH_min_vaf = 0.05                           # lower limit for VAF included in MATH score
MATH_max_vaf = 0.75                           # upper limit for VAF included in MATH score



####################
## Nothing below this should need editing
####################

# Define output directory for filtered VCF files

## Extract list of files to analyze
filenames <- list.files(path=path_to_files, pattern = "*.vcf.gz$", recursive = recursive_value, ignore.case = FALSE, full.names=TRUE)
filenames
found <- grepl(pattern_to_match, filenames, ignore.case = TRUE)  # find e.g. normal or tumor files only
filenames <- filenames[found]
filenames
filenames_all <- filenames
basenames <- basename(filenames)
basenames

# Read a file of tumor stage or other meta data
samples.data <- read.csv(samples_file, header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
samples.data

# Keep only files marked include = yes
if (include_marked_only == TRUE) {
    inc <- samples.data$INCLUDE == "yes"
    samples.data <- samples.data[inc,]
}
samples.data

keep <- samples.data$PRIMARY_FILE %in% basenames
summary(keep)
samples.data <- samples.data[keep,]

keep <- basenames %in% samples.data$PRIMARY_FILE
summary(keep)

basenames <- basenames[keep]
filenames <- filenames[keep]

ff <- t(data.frame(sapply(basenames, strsplit, split='_')))
ff
samplenames <- as.matrix(paste(ff[,1], ff[,2], sep="_"))
samplenames

# set up data file
ms_data <- data.frame(SAMPLES=samplenames, MATH=0, VARIANTS=0)

basenames
filenames

no.samples <- length(samplenames)
no.samples

# make sure the order of files and meta data is correct
m <- match(basenames, samples.data$PRIMARY_FILE)
m
samples.data <- samples.data[m,]

labels <- samplenames
labels

# Load all the libraries and functions
source("R/math_score_functions.R", echo=TRUE)

outputdir = paste(path_to_save, "/", analysis.name, "_", analysis.date, "/", sep="")
outputdir

mkdirs(outputdir)                           # Creates output directory if it doesn't exist, old data gets over-written if it does
testDirectory(outputdir, "rw")              # Test if output directory exists and is writeable
stopifnot(testDirectory(outputdir, "rw"))   # stop script if output directory not writeable
outputdir


for (s in 1:no.samples) {
    source("R/math_score_script.R", echo=TRUE)
}

# Save ms data file
outputfile = paste(outputdir, "ms_data_stage.csv", sep="")
write.table(ms_data, file=outputfile, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE, qmethod = "double")

# Add ms data to samples.data file
samples.data$MATH <- ms_data$MATH
samples.data$VARIANTS <- ms_data$VARIANTS

# Save complete sample data file
outputfile = paste(outputdir, "step1_sample_data.csv", sep="")
write.table(samples.data, file=outputfile, quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE, qmethod = "double")




####################
## Prepare Stage II MATH score boxplot
####################

samples.data

# Perform t.test of different groups
group2 <- samples.data$GROUP_NO == 1
group2 <- samples.data[group2,]
group2
no.group2 <- length(group2[,1])

group3 <- samples.data$GROUP_NO == 2
group3 <- samples.data[group3,]
group3
no.group3 <- length(group3[,1])


ms.stage.2v3.t.test <- t.test(as.numeric(as.matrix(group2$MATH)), as.numeric(as.matrix(group3$MATH)), var.equal = TRUE)
ms.stage.2v3.t.test


plot(samples.data$GROUP_NO, (samples.data$MATH))


# Make a boxplot
y = as.numeric(as.matrix(samples.data$MATH))
grp = samples.data$GROUP_NO
title = "Tumor Heterogeneity (MATH) in Stage II Tumors"
xlabel = "Metastasis Status"
ylabel = "MATH Score (100*MAD/median VAF at heterozygous tumor loci)"
xlab1 <- paste("No Metastasis (n=", no.group2, ")", sep="")
xlab2 <- paste("Metastasis Positive (n=", no.group3, ")", sep="")
grp_names = c(xlab1, xlab2)
p_val1 <- paste("pval=", round(ms.stage.2v3.t.test$p.val, digits = 3), sep="")
y_min <- round(min(samples.data$MATH)-1,0)
y_max <- round(max(samples.data$MATH)+3, 0)
ylimits <- c(y_min,y_max)
borders = c("#000000", "#FF0000")
colors = c("#FFFFFF", "#ffffff")


name = paste(outputdir, "MATH_score_stageII_boxplot.jpg",  sep = "", collapse = NULL)
jpeg(name, height = 700, width = 700)
boxplot(y ~ grp, main = title, xlab = xlabel, ylab = ylabel, ylim = ylimits, varwidth = TRUE, names = grp_names, border = borders, col = colors)
par(xpd=TRUE)
yrange<-par("usr")[3:4]
ypos1<-(yrange[2]+diff(yrange)/40) - 3
segments(1,ypos1,2,ypos1)
text(1.5,(1+ypos1+diff(yrange)/40),p_val1)
if (ms.stage.2v3.t.test$p.val < 0.05) {
    text(1.5,(ypos1+diff(yrange)/40),"*",cex=2)
}
par(xpd=FALSE)
dev.off()

name = paste(outputdir, "MATH_score_stageII_boxplot.pdf",  sep = "", collapse = NULL)
pdf(name, height = 7.5, width = 7.5)
boxplot(y ~ grp, main = title, xlab = xlabel, ylab = ylabel, ylim = ylimits, varwidth = TRUE, names = grp_names, border = borders, col = colors)
par(xpd=TRUE)
yrange<-par("usr")[3:4]
ypos1<-(yrange[2]+diff(yrange)/40) - 3
segments(1,ypos1,2,ypos1)
text(1.5,(1+ypos1+diff(yrange)/40),p_val1)
if (ms.stage.2v3.t.test$p.val < 0.05) {
    text(1.5,(ypos1+diff(yrange)/40),"*",cex=2)
}
par(xpd=FALSE)
dev.off()


# Save complete session info
sessionfile = paste(outputdir, "Session_Info.txt", sep="")
si <- capture.output(sessionInfo())
write.table(si, file = sessionfile, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)


