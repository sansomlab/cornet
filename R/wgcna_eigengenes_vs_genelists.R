stopifnot(
  require(optparse),
  require(openxlsx),
  require(WGCNA),
  require(ggplot2),
  require(dendsort),
  require(RColorBrewer),
  require(ComplexHeatmap)
)

# Options ----
rundir <- "/gfs/work/ssansom/covid_spatial/wgcna_covid_qn"
option_list <- list(
  make_option(
    c("--input"),
    default=file.path(rundir,
                      "qn.wgcna.dir/clean.dir/clean.RData"),
    help="where should the output files be saved"
  ),
  make_option(
    c("--outdir"),
    default="test/eg_vs_genelists.dir",
    help="where should the output files be saved"
  ),
  make_option(
    c("--annotation"),
    default="test/data.dir/anno.data.tsv",
    help='A file containing "gene_id" and "gene_name" columns'
  ),
  make_option(
    c("--modules"),
    default=file.path(rundir,
                      "qn.wgcna.dir/modules.dir/modules.RData"),
    help='the input data, genes in rows, samples in columns. A column of identifers, named "gene_id" is required.'
  ),
  make_option(
    c("--genelists"),
    default="covid.spatial.genelists.tsv",
    help='the genelists with columns gene_name, category, plot_group'
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
    c("--threads"),
    default=2,
    help='number of WGCNA threads'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))
message("running with options")
print(opt)

# make the output directory if it does not exist
dir.create(opt$outdir,
           showWarnings = FALSE, recursive=TRUE)

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


clustFunc <- function(x)
{
  x <- as.dist(1-cor(as.matrix(x)))
  x <- dendsort(hclust(x))
  x
}

genelistHeatmap <- function(datExpr,
                            genelists,
                            annotation=anno,
                            idcol="ensembl_id",
                            namecol="gene_name",
                            plot_group=NULL)
{
  if(!is.null(plot_group))
  {
    genelist <- genelists[genelists$plot_group==plot_group,]
  } else { genelist <- genelists }

  categories <- unique(genelist$category)

  if(length(categories)>1)
  {

  cat_cols <- brewer.pal(length(categories),
                         "Dark2")[1:length(categories)]
  names(cat_cols) <- categories

  cann <- HeatmapAnnotation(category=genelist$category,
                            col=list("category"=cat_cols))

  } else { cann <- NULL }


  genelist <- genelist[genelist[[idcol]] %in% colnames(datExpr),]

  x <- datExpr[, genelist[[idcol]]]
  colnames(x) <- annotation[colnames(x), namecol]

  moduleTraitCor = cor(MEs, x, use = "p")

  Heatmap(moduleTraitCor,
          top_annotation = cann,
          cluster_rows = clustFunc(t(moduleTraitCor)),
          cluster_columns = clustFunc(moduleTraitCor),
          column_names_gp = grid::gpar(fontsize = 8))
}

genelists <- read.table(opt$genelists,
                        header=T, sep="\t", as.is=T)

# subset to genes present in the dataset
genelists <- genelists[genelists[[opt$idcol]] %in% colnames(datExpr),]

for(plot_group in unique(genelists$plot_group))
{
  print(plot_group)
  png(file.path(opt$outdir,
                paste("genelist",plot_group,"png", sep=".")),
                width=10, height=6, unit="in", res=300)

  x <- genelistHeatmap(datExpr,
                       genelists,
                       annotation=anno,
                       idcol=opt$idcol,
                       namecol=opt$namecol,
                       plot_group=plot_group)
  print(x)
  dev.off()
}
