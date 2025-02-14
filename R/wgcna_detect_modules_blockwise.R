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
  require(fastcluster),
  require(ggplot2)
)

# Options ----

option_list <- list(
  make_option(
    c("--input"),
    default="test/clean.dir/dataInput.RData",
    help='The clean data in an RData file.'
  ),
  make_option(
    c("--outdir"),
    default="test/modules.dir",
    help="where should the output files be saved"
  ),
  make_option(
    c("--outfilename"),
    default="modules.RData",
    help="The name for the file containing the processed data"
  ),
  make_option(
    c("--threads"), default=4,
    help='Number of threads for parallel operations.'
  ),
  make_option(
    c("--maxblocksize"), default=Inf,
    help='Number of blocks to use'
  ),
  make_option(
    c("--softpower"),
    default=4,
    help="The soft thresholding power"
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
  ),
  make_option(
    c("--tomtype"),
    default="unsigned",
    help="The TOM Type",
  ),
  make_option(
    c("--minmodulesize"),
    default=30,
    help="minimum number of genes in a module"
  ),
  make_option(
    c("--medissthreshold"),
    default=0.25,
    help='dissimilarity threshold for merging modules'
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


# -------------------- 2. Adjacency -------------------------- #

if(opt$adjcorfnc=="pearson")
{
  corfnc = "cor"
  coropt = list(use = "p")
} else if(opt$adjcorfnc=="spearman")
{
  corfnc = "cor"
  coropt = list(use = "p", method = "spearman")
}


bwnet = blockwiseModules(datExpr,
                         maxBlockSize = opt$maxblocksize,
                         power = opt$softpower,
                         corType = opt$adjcorfnc,
                         networkType = opt$networktype,
                         TOMType = opt$tomtype,
                         minModuleSize = opt$minmodulesize,
                         reassignThreshold = 0,
                         mergeCutHeight = opt$medissthreshold,
                         numericLabels = TRUE,
                         saveTOMs = FALSE,
#                         saveTOMFileBase = "femaleMouseTOM-blockwise",
                         verbose = 3)


# Plot the block dendrograms.
nblocks = length(bwnet$dendrograms)

# The merged module colors
moduleLabels = bwnet$colors;
moduleColors = labels2colors(moduleLabels)
# Eigengenes of the merged modules:
MEs = bwnet$MEs;

me_col <- labels2colors(as.numeric(gsub("ME","",colnames(MEs))))
colnames(MEs) <- paste0("ME",me_col)

# visualise the blocks
for(block in 1:nblocks)
{
  sizeGrWindow(6,6)
  pdf(file.path(opt$outdir, paste0("block ",block, " module_dendrogram.pdf")),
      width = 6, height = 6)
  plotDendroAndColors(bwnet$dendrograms[[block]], moduleColors[bwnet$blockGenes[[block]]],
                      "Module colors", main = paste0("Gene dendrogram and module colors in block", block),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}


# ---------- 5. Visualise the similarity of the modules ----------- #

sizeGrWindow(6,6)

pdf(file.path(opt$outdir, "eigengene_dendrogram.pdf"),
    width=6, height=6)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()

pdf(file.path(opt$outdir, "eigengene_adjacency_heatmap.pdf"),
    width=6, height=6)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

# ------------------ 6. Fix colors and save ----------------------- #

# Here we will save only
# 1. the module labels
# 2. the module eigen gene matrix (MEs)

save(MEs, moduleColors,              # moduleColors, geneTree, moduleLabels,
     file = file.path(opt$outdir,
                      opt$outfilename))
