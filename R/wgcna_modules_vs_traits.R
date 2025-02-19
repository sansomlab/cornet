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

options(bitmapType = "cairo")

# Options ----

option_list <- list(
  make_option(
    c("--input"),
    default="test/clean.dir/dataInput.RData",
    help='An RData file containing the clean data.'
  ),
  make_option(
    c("--modules"),
    default="test/modules.dir/modules.RData",
    help='the RData file containing the MEs, moduleLabels and moduleColors objects.'
  ),
  make_option(
    c("--annotation"),
    default="test/data.dir/anno.data.tsv",
    help='A file containing "gene_id" and "gene_name" columns'
  ),
  make_option(
    c("--idcol"),
    default="gene_id",
    help='the column containing the gene identifiers'
  ),
  make_option(
    c("--namecol"),
    default="gene_name",
    help='the column containing the gene names'
  ),
  make_option(
    c("--outdir"),
    default="test/traits.dir",
    help="where should the output files be saved"
  ),
  make_option(
    c("--outfilename"),
    default="gene_info.tsv",
    help="The name for the file containing the processed data"
  ),
  make_option(
    c("--threads"), default=4,
    help='Number of threads for parallel operations.'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))
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

# Load the module assignment information
mnames = load(file = opt$modules);
#The variable lnames contains the names of loaded variables.
message("Names of the loaded module data objects: ")
print(mnames)

# Load the annotation
anno = read.table(gzfile(opt$annotation), header=T, as.is=T,
                   sep="\t")


anno <- anno[, c(opt$idcol, opt$namecol)]
anno <- unique(anno)
# if ids map to multiple symbols, only the first given symbol will be used.
rownames(anno) <- make.unique(anno[[opt$idcol]])

# -------------------- 2. trait heatmap ----------------- #

# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels
# MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs)

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

pdf(file.path(opt$outdir, "module_trait_relationships.pdf"),
    width=10, height=6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
## Display the correlation values within a heatmap plot

message("names dat traits")
print(head(datTraits))
print(names(datTraits))
message("names me")
print(names(MEs))


labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# ------------------ 3. gene module membership -------------------------

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

# ------------------------- 4. write out a summary table ------------------------ #

ids <- names(datExpr)

mods <- unique(moduleColors)

# order modules from smallest to largest
modOrder <- order(table(moduleColors))

# Create the starting data frame
geneInfo0 = data.frame(gene_id = anno[ids,opt$idcol],
                       gene_name = anno[ids, opt$namecol],
                       moduleColor = moduleColors)
#                       geneTraitSignificance,
#                       GSPvalue)


# geneModuleMembership
MMPvalue <- MMPvalue[rownames(geneModuleMembership),]

# Order modules by their significance for weight
# modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
begin <- TRUE
for (xname in mods[modOrder])
{
  genes_in_module = colnames(datExpr)[moduleColors==xname]

  data <- data.frame(membership=geneModuleMembership[genes_in_module, paste("MM", xname, sep="")],
                     p.value=MMPvalue[genes_in_module, paste("p.MM",xname,sep="")],
                     gene_id=genes_in_module,
                     gene_name=anno[genes_in_module, opt$namecol])
  data$module <- xname

  data <- data[order(data$p.value),]
  data <- data[, c("module", "gene_id", "gene_name", "membership", "p.value")]

  if(begin) { gene_info <- data
  begin <- FALSE
  } else { gene_info <- rbind(gene_info, data) }


}


write.table(gene_info, file.path(opt$outdir, opt$outfilename),
            sep="\t", col.names=T, row.names=F, quote=F)


# Write out a table containing the eigengenes.
xx <- MEs
xx_names <- colnames(xx)
xx$sample_id <- rownames(xx)
xx <- xx[,c("sample_id",xx_names)]

write.table(xx, file.path(opt$outdir, "eigengenes.tsv"),
            sep="\t", col.names=T, row.names=F, quote=F)
