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
    help='The clean data RData file.'
  ),
  make_option(
    c("--outdir"),
    default="test/modules.dir",
    help="where should the output files be saved"
  ),
  make_option(
    c("--outfilename"),
    default="adjacency.RData",
    help="The outfile name"
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


# -------------------- 2. Adjacency -------------------------- #

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


# >>>>>>>>>>>>>>>>>>>>>>>>>>> stepwise
# >>>>>>>>>>>>>>>>>>>>>>>>>>> or blockwise


adjMat <- adjacency(datExpr,
                    type = opt$networktype,
                    power = opt$softpower,
                    corFnc = corfnc,
                    corOptions = coropt,
                    distFnc = "dist")
#                      disopt$adjdistfnc)

dim(adjMat)


# ------------------ 3. asve ----------------------- #

# Here we will save only
# 1. the module labels
# 2. the module eigen gene matrix (MEs)

save(adjMat,              # moduleColors, geneTree, moduleLabels,
     file = file.path(opt$outdir,
                      opt$outfilename))
