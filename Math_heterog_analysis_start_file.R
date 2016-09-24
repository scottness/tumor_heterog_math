# Set of R scripts for analyzing math scores
# Scott Ness, Sep 2016 sness@salud.unm.edu
# 
#### Description
# 
# Important scripts are in the 'R' folder
# 

# Update the information in the next few lines, then run...
#
## Define the file locations and analysis parameters
analysis.name = "math_analysis"      # Define this for each project or set of samples
path_to_save = "/Users/scottness/Dropbox/In\ Preparation/Heterogeneity_paper/data_analysis"   

# Location, files to be analyzed and parameters
# path_to_files = "/Volumes/Genomics1/data_lake/Colon_paper_2016/ccp_orig_vcf_tsvc4421"   # work
path_to_files = "/Volumes/Data1/data_lake/Colon_paper_2016/ccp_orig_vcf_tsvc4421"  # home
vcfs_from_GATK = FALSE              # TRUE if analyzing RNA-seq VCF files
recursive_value = FALSE              # TRUE to look through sub-folders
pattern_to_match = ".fixed.vcf.gz$"   # All the files must have this in their name

# samples_file = "/Volumes/Genomics1/data_lake/Colon_paper_2016/colon_TvN_sample_data.csv"   # work
samples_file = "/Volumes/Data1/data_lake/Colon_paper_2016/colon_math_sample_data.csv"   # home
# Samples data should be a tab-delimited text file with these fields:
# UID               (file ID, part of sample name)
# INCLUDE           (yes or no)
# GROUP_LETTER      (Letter code: A,B,C,D,E)
# GROUP_NO          (Number: 1,2,3)
# PRIMARY_FILE      (tumor file name)
# Additional columns are allowed

include_marked_only = TRUE        # If true, only use samples with "yes" in the samples data INCLUDE column

# s=0

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


# Load gene lists and database locations (e.g. COSMIC)
source("../lists_locations_hg19.R")


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
## Compare Additional Stage II samples
####################
rm(list=ls())


# Update the information in the next few lines, then run...
#
## Define the file locations and analysis parameters
analysis.name = "math_analysis_stageII"      # Define this for each project or set of samples
path_to_save = "/Users/scottness/Dropbox/In\ Preparation/Heterogeneity_paper/data_analysis"   

# Location, files to be analyzed and parameters
# path_to_files = "/Volumes/Genomics1/data_lake/Colon_paper_2016/ccp_orig_vcf_tsvc4421"   # work
path_to_files = "/Volumes/Data1/data_lake/Colon_paper_2016/ccp_stgII_vcf_tsvc"  # home
vcfs_from_GATK = FALSE              # TRUE if analyzing RNA-seq VCF files
recursive_value = FALSE              # TRUE to look through sub-folders
pattern_to_match = ".fixed.vcf.gz$"   # All the files must have this in their name

# samples_file = "/Volumes/Genomics1/data_lake/Colon_paper_2016/colon_TvN_sample_data.csv"   # work
samples_file = "/Volumes/Data1/data_lake/Colon_paper_2016/colon_math_stageII_sample_data.csv"   # home
# Samples data should be a tab-delimited text file with these fields:
# UID               (file ID, part of sample name)
# INCLUDE           (yes or no)
# GROUP_LETTER      (Letter code: A,B,C,D,E)
# GROUP_NO          (Number: 1,2,3)
# PRIMARY_FILE      (tumor file name)
# Additional columns are allowed

include_marked_only = TRUE        # If true, only use samples with "yes" in the samples data INCLUDE column

# s=0

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


# Load gene lists and database locations (e.g. COSMIC)
source("../lists_locations_hg19.R")


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
outputfile = paste(outputdir, "ms_data.csv", sep="")
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


