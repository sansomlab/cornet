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
  require(gsfisher),
  require(ggplot2)
)

# Options ----

option_list <- list(
  make_option(c("--input"),
              default="test/traits.dir/membership.tsv",
              help='Table of module gene memberhsip. A column of identifers, named "gene_id" is required.'),
  make_option(c("--outdir"), default="test.dir/genesets.dir",
              help="outdir"),
  make_option(c("--module"),
              default="turquoise",
               help='Module to run geneset enrichment on'),
  make_option(c("--species"), default="hs",
              help="species: hs or mm"),
  make_option(c("--annotation"),default="none",
              help="entrez_id,ensembl_id,gene_name tab sep"),
  make_option(c("--idcol"), default="gene_id",
              help="the column of the annotations file containing the gene ids"),
  make_option(c("--kegg_pathways"), default="none",
              help="R object containing kegg_pathways list"),
  make_option(c("--gmt_files"), default="none",
              help="comma separated list of gmt files"),
  make_option(c("--gmt_names"), default="none",
              help="comma separated list of names for the gmt files"),
  make_option(c("--project"), default="WgcnaAnalysis",
              help="project name"),
  make_option(c("--prefix"), default="genesets",
              help="prefix for out files")

)

opt <- parse_args(OptionParser(option_list=option_list))
message("Running with options:")
print(opt)


# make the output directory if it does not exist
dir.create(opt$outdir, 
           showWarnings = FALSE)

# --------------------- 1. setup ---------------------- #

message("retrieving the gene lists")

mm <- read.table(opt$input, sep="\t", header=T,
                         as.is=T)

foreground = mm$gene_id[mm$module == opt$module]
universe = unique(mm$gene_id)

## set the run specs
run_specs <- paste(opt$numpcs,opt$resolution,opt$algorithm,opt$testuse,sep="_")
#background <- read.table(gzfile(opt$universe), header=T, as.is=T, sep="\t")$gene_id

if(length(universe)<100)
{
  stop("Less than 100 genes in the gene universe!!")
}

print(paste0("no ensembl_ids in foreground: ", length(foreground)))
print(paste0("no ensembl_ids in universe: ", length(universe)))

message("loading annotations")
anno <- read.table(gzfile(opt$annotation), as.is=T, sep="\t", header=T)
print(head(anno))

## get the foreground and universe gene lists
fg_entrez <- unique(anno$entrez_id[anno[[opt$idcol]] %in% foreground])
u_entrez <- unique(anno$entrez_id[anno[[opt$idcol]] %in% universe])

fg_entrez <- as.character(fg_entrez[!is.na(fg_entrez)])
u_entrez <- as.character(u_entrez[!is.na(u_entrez)])

print(paste0("no entrez_ids in foreground: ", length(fg_entrez)))
print(paste0("no entrez_ids in universe: ", length(u_entrez)))

message("performing enrichment tests")
## perform the enrichment tests
if(length(fg_entrez)>0)
{
  outPrefix = paste0(opt$outdir,"/", opt$prefix, ".",opt$module)
  kegg_pathways <- readRDS(opt$kegg_pathways)
  
  gmts <- list()
  
  if(!opt$gmt_names == "none" | !opt$gmt_files == "none")
  {
    ## parse the list of gmts
    gmt_names <- strsplit(opt$gmt_names,",")[[1]]
    gmt_files <- strsplit(opt$gmt_files,",")[[1]]
    
    if(length(gmt_names) != length(gmt_files))
    {
      stop("the same number of gmt_names and gmt_files must be specified")
    }
    
    for(i in 1:length(gmt_names))
    {
      gmts[[gmt_names[i]]] <- gmt_files[i]
    }
  }
  
  ## run the analysis.
  results <- analyseGenesets(fg_entrez,
                             u_entrez,
                             gene_id_type="entrez",
                             kegg_pathways=kegg_pathways,
                             gmt_files=gmts,
                             species=opt$species)
  
  ## outPrefix,
  for(geneset in names(results))
  {
    write.table(results[[geneset]],
                gzfile(paste(outPrefix,geneset,"tsv","gz",sep=".")),
                row.names=FALSE, col.names=TRUE, quote=FALSE,
                sep="\t")
  }
  
} else {
  print("module has no significant genes")
}

message("enrichment tests complete")
