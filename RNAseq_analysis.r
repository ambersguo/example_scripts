library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)


# Read the data into R
seqdata <- read.delim("RawCountFile_rsemgenes_for_RNAseq_analysis.txt", stringsAsFactors = FALSE)
# Read the sample information into R
sampleinfo <- read.delim("sample_info.txt")
# Remove the first column from seqdata
countdata <- seqdata[,-1]
# Look at the output
head(countdata)

# Store gene id as rownames
rownames(countdata) <- seqdata[,1]
# Look at the output
head(countdata)

# Obtain CPMs
myCPM <- cpm(countdata)
# Have a look at the output
head(myCPM)

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
table(rowSums(thresh))

# we would like to keep genes that have 2 TRUES in each row of thresh
keep <- rowSums(thresh) == 2
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
summary(keep)
# check how many genes left
dim(counts.keep)

# convert counts to DGEList object
y <- DGEList(counts.keep)
y
# See what slots are stored in y
names(y)
# Library size information is stored in the samples slot
y$samples

pdf("barplot_of_library_sizes.pdf")
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(y$samples$lib.size,names=colnames(y),las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()

pdf("boxplot_of_logCPMs.pdf")
# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalized)")
dev.off()


########## Hierarchical clustering with heatmaps #########

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var_100 <- names(sort(var_genes, decreasing=TRUE))[1:100]
select_var_1000 <- names(sort(var_genes, decreasing=TRUE))[1:1000]

# Subset logcounts matrix
highly_variable_lcpm_100 <- logcounts[select_var_100,]
highly_variable_lcpm_1000 <- logcounts[select_var_1000,]

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
         
# Plot the heatmap
pdf(file="Heatmap_top100_high_var_genes_with_labels.pdf")
heatmap.2(highly_variable_lcpm_100,col=rev(morecols(50)),trace="none",
          scale="none", symm=F,symkey=F,symbreaks=F,
         margins=c(8,8), cexRow=0.2, cexCol=1)
dev.off()

pdf(file="Heatmap_top1000_high_var_genes.pdf")
heatmap.2(highly_variable_lcpm_1000,col=rev(morecols(50)),trace="none",
          scale="none", symm=F,symkey=F,symbreaks=F, labRow=F,
          margins=c(8,8), cexRow=0.2, cexCol=1)
dev.off()
