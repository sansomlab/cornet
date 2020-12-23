## Title ----
##
## Make plots visualising input data so that it can be cleaned.
##
## Description ----
##
## Follows this tutorial: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf
##
## Details ----
##
##
## Usage ----
##
## See options.

# Libraries ----

stopifnot(
  require(WGCNA),
  require(ggplot2),
  require(optparse),
  require(fastcluster)
)

# Options ----

option_list <- list(
  make_option(
    c("--outdir"),
    default="test/clean.dir",
    help="where should the output files be saved"
  ),
    make_option(
    c("--input"),
    default="test/data.dir/exprs.data.tsv",
    help='the input data, genes in rows, samples in columns. A column of identifers, named "gene_id" is required.'
    ),
  make_option(
    c("--idcol"),
    default="gene_id",
    help='the column containing the gene identifiers'
  ),
  make_option(
    c("--traitdata"),
    default="test/data.dir/trait.data.tsv",
    help='Table containing the trait data. Must contain column "sample_name"'
  ),
  make_option(
    c("--outfilename"),
    default="dataInput.RData",
    help="The name for the file containing the processed data"
  ),
  make_option(
    c("--minfraction"), default=0.5,
    help='minimum fraction of non-missing samples for "good" genes'
  ),
  make_option(
    c("--minnsamples"),
    default=4,
    help="The minimum number of non-missing samples for gene to be considered good."
  ),
  make_option(
    c("--minngenes"),
    default=4,
    help="miminum number of good genes for dataset to be considered good."
  ),
  make_option(
    c("--minrelativeweight"),
    default=0.1,
    help="observations whose relative weight is below this threshold are considered missing (see WGCNA docs)"
  ),
  make_option(
    c("--cutheight"),
    default=15,
    help="The cut height for excluding samples based on the sample dendrogram"
  ),
  make_option(
    c("--minsize"),
    default=50,
    help="minimum number of object on a branch to be considered a cluster"
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

message("Running with options:")
print(opt)


# make the output directory if it does not exist
dir.create(opt$outdir,
           showWarnings = FALSE)

# --------------------- 1. prepare input data ---------------------- #

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#Read in the input data set
input_data = read.table(opt$input, sep="\t",
                        header=T, as.is=T);

# Set the row names, ensuring that they are unique
rownames(input_data) <- make.unique(input_data[[opt$idcol]])
input_data[[opt$idcol]] <- NULL

# Take a quick look at what is in the data set:
message("Dimensions of input data: ", dim(input_data))
message("Names of input data: ", names(input_data))


# transpose the data
exprs_data = as.data.frame(t(input_data))


# ------------------------- 2. get good samples and genes --------------------- #

gsg = goodSamplesGenes(exprs_data,
                       weights = NULL,
                       minFraction = opt$minfraction,
                       minNSamples = opt$minnsamples,
                       minNGenes = opt$minngenes,
                       tol = NULL,
                       minRelativeWeight = opt$minrelativeweight,
                       verbose = 3);

message("All genes have passed the cuts? ", gsg$allOK)

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(exprs_data)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(exprs_data)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  exprs_data = exprs_data[gsg$goodSamples, gsg$goodGenes]
}

# ------------------- 3. visualise the sample tree --------------------- #


sampleTree = hclust(dist(exprs_data), method = "average");

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.


sizeGrWindow(12,9)

pdf(file = file.path(opt$outdir,
                     "sampleClustering.pdf"),
    width = 12, height = 9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = opt$cutheight, col = "red")


dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree,
                     cutHeight = opt$cutheight,
                     minSize = opt$minsize)

message("number samples above and below the cut:")

xx <- table(clust)

print(xx)

if(sum(as.numeric(names(xx)))==0)
{
  stop("No samples made the cut and/or number of samples in cluster (minsize) too small")
}



# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = exprs_data[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


# ---------------- 4. visualise the trait data ----------------------- #

if(!is.null(opt$traitdata))
{

traitData = read.table(opt$traitdata,
                       header=T,
                       sep="\t")

rownames(traitData) <- traitData$sample_name

traitData$sample_name<-NULL

message("Dimensions of the trait data: ")
print(dim(traitData))
message("Names of the trait data: ", names(traitData))

# Form a data frame analogous to expression data that will hold the clinical traits.
datTraits = traitData[rownames(datExpr),]

collectGarbage();


pdf(file = file.path(opt$outdir,
                     "sampleClustering_with_traits.pdf"),
    width = 12, height = 9)

par(cex = 0.6)
par(mar = c(0,4,2,0))

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

dev.off()

save(datExpr, datTraits, file = file.path(opt$outdir,
                                          opt$outfilename))

} else {

  save(datExpr, file = file.path(opt$outdir,
                                 opt$outfilename))
}
