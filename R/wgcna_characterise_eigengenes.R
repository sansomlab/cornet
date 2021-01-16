stopifnot(
  require(WGCNA),
  require(ggplot2),
  require(optparse),
  require(ComplexHeatmap),
  require(ggplot2),
  require(openxlsx),
  require(circlize),
  require(dendsort),
  require(yaml),
  require(dplyr)
)

options(bitmapType = "cairo")

# Options ----
rundir <- "/gfs/work/ssansom/covid_spatial/covid_old/wgcna_covid_qn"
option_list <- list(
  make_option(
    c("--outdir"),
    default="test/clean.dir",
    help="where should the output files be saved"
  ),
  make_option(
    c("--eigengenes"),
    default=file.path(rundir,
                      "qn.wgcna.dir/membership.dir/eigengenes.tsv"),
    help='the input data, genes in rows, samples in columns. A column of identifers, named "gene_id" is required.'
  ),
  make_option(
    c("--membership"),
    default=file.path(rundir,
                      "qn.wgcna.dir/membership.dir/membership.tsv"),
    help='the input data, genes in rows, samples in columns. A column of identifers, named "gene_id" is required.'
  ),
  make_option(
    c("--params"),
    default=file.path(rundir,
                      "pipeline.yml"),
    help='the pipeline yml file'
  ),
  make_option(
    c("--namecol"),
    default="gene_name",
    help='the column containing the gene names'
  ),
  make_option(
    c("--traitdata"),
    default=NULL,
    help='Table containing the trait data. Must contain column "sample_name"'
  ),
  make_option(
    c("--metadata"),
    default=NULL
    # default=file.path(rundir,"data.dir/qn.dir/meta.data.tsv"),
    # help='Table containing the trait data. Must contain column "sample_name"'
  ),
  make_option(
      c("--figwidth"),
      default=10,
      help="figure width in inches"),
  make_option(
      c("--figheight"),
      default=8,
      help="figure width in inches")

)

opt <- parse_args(OptionParser(option_list=option_list))

eg <- read.table(opt$eigengenes,
                 header=T, sep="\t", as.is=T)
rownames(eg) <- eg$sample_id
eg$sample_id <- NULL

if(!is.null(opt$traitdata))
{
trait_data <- read.table(opt$traitdata,
                         header=T, sep="\t", as.is=T)

rownames(trait_data) <- trait_data$sample_name
}

if(!is.null(opt$metadata))
{
  meta_data <- read.table(opt$metadata,
                          header=T, sep="\t", as.is=T)

  rownames(meta_data) <- meta_data$sample_name

  if(!is.null(opt$traitdata))
  {
      trait_data <- cbind(trait_data, meta_data[rownames(trait_data),])
  } else { trait_data <- meta_data }

}

params <- read_yaml(opt$params)

## 1. make the trait/metadata-annotated eigengene heatmap.

cfuns <- list()
for(trait in names(params$traits))
{
  message("setting color palette for trait: ", trait)
  if(is.list(params$traits[[trait]]))
  {
    cfuns[[trait]] <- unlist(params$traits[[trait]])
  } else {
  color_str <- params$traits[[trait]]
  color_str <- gsub(" ","", color_str)
  colors <- strsplit(color_str,",")[[1]]
  print(colors)
  ncol <- length(colors)
  quantiles <- c(0,1:(ncol-2)/(ncol-1),1)
  print(quantiles)
  breaks <- quantile(trait_data[[trait]],quantiles, na.rm=TRUE)
  print(breaks)
  cfuns[[trait]] <- colorRamp2(breaks,colors)
  }
}

trait_df <- trait_data[rownames(eg),names(params$traits)]

cann <- HeatmapAnnotation(df=trait_df,
                          col=cfuns)

corClust <- function(x)
{
  x <- as.dist(1-cor(as.matrix(x)))
  x <- dendsort(hclust(x))
  x
}

eg_mat <- t(as.matrix(eg))

## remove the grey module which contains the
## unassigned genes.
eg_mat <- eg_mat[!rownames(eg_mat)=="MEgrey",]

hm <- Heatmap(eg_mat,
        cluster_rows = corClust(t(eg_mat)),
        cluster_columns = corClust(eg_mat),
        top_annotation = cann)

png(file.path(opt$outdir,"eigengene_heatmap.png"),
    width=opt$figwidth,
    height=opt$figheight,
    units="in",
    res=300)
print(hm)
dev.off()

save(eg, eg_mat, trait_df, cann, corClust, hm, cfuns,
     file=file.path(opt$outdir, "eigengene.heatmap.Rdata"))


## 2. make the summary plot of the eigengene correlations.

membership <- read.table(opt$membership, header=T,
                         sep="\t", as.is=T)

membership <- membership[membership$membership > 0,]

modsizes <- table(membership$module)
modsizes <- modsizes[rev(order(modsizes))]

modsize_df <- data.frame(module=names(modsizes),
                         n=as.vector(modsizes))

modsize_df$module <- factor(modsize_df$module, levels=modsize_df$module)

cvals <- as.vector(modsize_df$module)
names(cvals) <- as.vector(modsize_df$module)

gp <- ggplot(modsize_df, aes(module, n, fill=module))
gp <- gp + scale_fill_manual(values=cvals)
gp <- gp + geom_bar(stat="identity", alpha=0.7, color="black")
gp <- gp + theme_classic() + ylab("numbers of genes positively\ncorrelated with module")
gp <- gp + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 legend.position = "none")

ggsave(file.path(opt$outdir,"eigengene_numbers.png"),
       width=6,height=3, units = "in")

membership$module <- factor(membership$module, levels=modsize_df$module)

gp <- ggplot(membership, aes(membership, y=(..count..), fill=module))
gp <- gp + scale_fill_manual(values=cvals)
gp <- gp + geom_density(alpha=0.7, color="black")
gp <- gp + facet_wrap(~module, ncol=6)
gp <- gp + theme_classic() + ylab("density of gene-eigengene correlations")
gp <- gp + xlab("correlation")
gp <- gp + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 strip.background = element_blank(),
                 legend.position= "none")

ggsave(file.path(opt$outdir,"eigengene_membership.png"),
       width=8,height=4.5, units = "in")

mdata <- membership %>%
  group_by(module) %>%
  top_n(n=20, wt=membership) %>%
  mutate(x = seq_along(membership))

mdata$y = 1
mdata$x = paste0("X",as.character(1:nrow(mdata)))
mdata$x <- factor(mdata$x, levels=mdata$x)

mdata <- data.frame(mdata)
mdata$module <- factor(mdata$module,levels=names(modsizes))

minm <- min(mdata$membership)
maxm <- max(mdata$membership)
breaks = minm + ((maxm - minm) * c(0,1/4,1/2,3/4,1))

gp <- ggplot(mdata, aes(x, y, fill=membership))
gp <- gp + geom_tile(stat="identity")
gp <- gp + scale_fill_gradientn(colors=c("grey","yellow","orange","brown"),
                                breaks=breaks,
                                labels=round(breaks,2),
                                name="correlation\nwith\neigengene")

gp <- gp + scale_x_discrete(breaks=mdata$x,
                            labels=mdata[[opt$namecol]])
gp <- gp + facet_wrap(~module, drop=T, ncol=3, scales="free_x")
gp <- gp + theme_classic()

gp <- gp + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 axis.title.y=element_blank(),
                 axis.title.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.text=element_text(size=8),
                 strip.background = element_blank(),
                 plot.margin = unit(c(0,0.5,0,0.8), "cm"),
                 panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(0, "lines"))

ggsave(file.path(opt$outdir,"eigengene_correlations.png"),
       width=10,height=10, units = "in")









# g <- ggplot_gtable(ggplot_build(gp))
# stripr <- which(grepl('strip-t', g$layout$name))
# fills <- rev(names(modsizes))
# k <- 1
# for (i in stripr) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# grid::grid.draw(g)
