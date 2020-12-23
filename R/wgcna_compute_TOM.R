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
    default="test/modules.dir/adjacency.RData",
    help='the R object containing the adjacency matrix'
  ),
  make_option(
    c("--outdir"),
    default="test/modules.dir",
    help="where should the output files be saved"
  ),
  make_option(
    c("--outfilename"),
    default="TOM.RData",
    help="The name for the file containing the processed data"
  ),
  make_option(
    c("--threads"), default=4,
    help='Number of threads for parallel operations.'
  ),
  make_option(
    c("--tomtype"),
    default="unsigned",
    help="The TOM Type",
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

opt$networktype <- gsub("_", " ", opt$networktype)

message("Running with options:")
print(opt)



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


# --------------------- 2. Topological overlap --------------------- #

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjMat,
                    TOMType = opt$tomtype)

dissTOM = 1-TOM


# ------------------ 3. Save the TOM ----------------------- #

# Here we will save only
# 1. the module labels
# 2. the module eigen gene matrix (MEs)

save(dissTOM,              # moduleColors, geneTree, moduleLabels,
     file = file.path(opt$outdir,
                      opt$outfilename))
