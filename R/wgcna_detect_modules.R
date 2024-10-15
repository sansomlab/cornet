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

options(bitmapType = "cairo")

# Options ----

option_list <- list(
  make_option(
    c("--cleandata"),
    default="test/clean.dir/dataInput.RData",
    help='the input data, genes in rows, samples in columns. A column of identifers, named "gene_id" is required.'
  ),
  make_option(
    c("--tomdata"),
    default="test/modules.dir/TOM.RData",
    help='the input data, genes in rows, samples in columns. A column of identifers, named "gene_id" is required.'
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
    c("--softpower"),
    default=4,
    help="The soft thresholding power"
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
  ),
  make_option(
    c("--deepsplit"),
    default=2,
    help='number from 0-4 for deep split parameter'
  ),
  make_option(
    c("--adjcorfnc"),
    default="_must_set_",
    help=paste("the function to be used to calculate co-expression similarity",
               "supported options are pearson, spearman")
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

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
lnames1 = load(file = opt$cleandata)
lnames2 = load(file = opt$tomdata)

## The variable lnames contains the names of loaded variables.
message("Names of the loaded variables: ")
print(lnames1)
print(lnames2)


# --------------------- 2. Cluster by topological overlap --------------------- #


# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

## Plot the resulting clustering tree (dendrogram)
pdf(file = file.path(opt$outdir,
                     "clustering_of_genes_by_topological_overlap.pdf"),
    width = 12, height = 9)

plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

dev.off()


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = opt$minmodulesize
dpsplit = opt$deepsplit


# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree,
                            distM = dissTOM,
                            deepSplit = dpsplit,
                            pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

message("table of modules")
print(table(dynamicMods))

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath

pdf(file = file.path(opt$outdir,
                     "module_dendrogram.pdf"),
    width = 12, height = 9)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

dev.off()

# -------------- 3. module merging ---------------------- #

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


# Calculate eigengenes
MEList = moduleEigengenes(datExpr,
                          colors = dynamicColors,
                          softPower = opt$softpower
                          )
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");


# Plot the result
pdf(file = file.path(opt$outdir,
                     "clustering_of_module_eigengenes.pdf"),
    width = 7, height = 6)

plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = opt$medissthreshold
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

dev.off()


# Call an automatic merging function
merge = mergeCloseModules(datExpr,
                          dynamicColors,
                          cutHeight = MEDissThres,
                          corFnc = corfnc,
                          corOptions = coropt,
                          verbose = 3)

# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
# moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs

pdf(file = file.path(opt$outdir,
                     "merged_module_dendrogram.pdf"),
    width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


# ---------- 4. Visualise the similarity of the modules ----------- #

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

# ------------------ 5. Fix colors and save ----------------------- #

# Here we will save only
# 1. the module labels
# 2. the module eigen gene matrix (MEs)

save(MEs, moduleColors,              # moduleColors, geneTree, moduleLabels,
     file = file.path(opt$outdir,
                      opt$outfilename))
