## Title ----
##
## Run step by step module detection
##
## Description ----
##
## With optional blockwise-approach for speed.
## 
## Details ----
##
##
## Usage ----
##
## See options.

# Libraries ----

stopifnot(
  require(optparse),
  require(WGCNA),
  require(ggplot2)
)

# Options ----

option_list <- list(
  make_option(
    c("--input"),
    default="test/clean.dir/dataInput.RData",
    help='the input data, genes in rows, samples in columns. A column of identifers, named "gene_id" is required.'
  ),
  make_option(
    c("--outdir"),
    default="test/soft.power.dir",
    help="where should the output files be saved"
  ),
  make_option(
    c("--threads"), default=4,
    help='Number of threads for parallel operations.'
  ),
  make_option(
    c("--networktype"),
    default="_must_set_",
    help="the type of network (unsigned, signed, signed hybrid or distance)"
  ),
  make_option(
    c("--adjcorfnc"),
    default="_must_set_",
    help=paste("the function to be used to calculate co-expression similarity",
               "supported options are pearson, spearman")
  ),
  make_option(
    c("--adjdistfnc"),
    default="_must_set_",
    help="distance function used when calculating co-expression similarity"
  )


)

opt <- parse_args(OptionParser(option_list=option_list))
opt$networktype <- gsub("_", " ", opt$networktype)

message("Running with options:")
print(opt)

if("_must_set_" %in% unlist(opt))
{
  stop("Critical (must set) parameter not specified")
}


# make the output directory if it does not exist
dir.create(opt$outdir, 
           showWarnings = FALSE)

# --------------------- 1. setup ---------------------- #

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads(nThreads=opt$threads)

# Load the data saved in the first part
lnames = load(file = opt$input);
#The variable lnames contains the names of loaded variables.
message("Names of the loaded variables: ")
print(lnames)


# -------------------- 2. soft-thresholding ----------------- #

if(opt$adjcorfnc=="pearson")
{
  corfnc = "cor"
  coropt = list(use = "p")
} else if(opt$adjcorfnc=="spearman")
{
  corfnc = "cor"
  coropt = list(use = "p", method = "spearman")
} else if(opt$adjcorfnc=="bicor")
{
  corfnc = "bicor"
  coropt = list(use = "p")
} else {
  
  stop("Correlation function not recognised")
}

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

#Call the network topology analysis function
sft = pickSoftThreshold(datExpr,
                        networkType = opt$networktype,
                        corFnc = corfnc,
                        corOptions = coropt,
                        powerVector = powers, 
                        verbose = 5)


# Plot the results:
sizeGrWindow(9, 5)

pdf(file = file.path(opt$outdir,"network_topology.pdf"),
    width = 9, height = 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()