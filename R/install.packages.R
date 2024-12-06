stopifnot(
    require(devtools),
    require(BiocManager))

install.package <- function(x, source="CRAN")
{
    if(source=="github")
    {
        xx <- strsplit(x,"/")[[1]]
        pkg <- xx[length(xx)]
        print(pkg)
        message(pkg)
    } else { pkg <- x }

    if (!require(pkg, character.only = TRUE))
    {
        message("installing missing package: ", pkg)

        if(source=="CRAN")
            {
                install.packages(pkg, dep=TRUE)

            } else if (source=="bioconductor")
            {
                BiocManager::install(pkg)

            } else if (source=="github")
            {
                devtools::install_github(x, upgrade="ask")
            }

        if(!require(pkg,character.only = TRUE)) stop("Package not found")

    }
}

cran_packages <- c("ape",
                   "circlize",
                   "colormap",
                   "dplyr",
                   "dendsort",
                   "ggplot2",
                   "gplots",
                   "fastcluster",
                   "openxlsx",
                   "optparse",
                   "RColorBrewer",
                   "reshape2",
                   "xtable",
                   "yaml")

bioconductor_packages <- c("ComplexHeatmap",
                           "WGCNA")

github_packages <- c("sansomlab/gsfisher")#,
#                     "sansomlab/tenx/tenxutils")

message("installing cran packages")
for(x in cran_packages)
{
    install.package(x, source="CRAN")
}

message("installing bioconductor packages")
for(x in bioconductor_packages)
{
    install.package(x, source="bioconductor")
}

message("installing github packages")
for(x in github_packages)
{
    install.package(x, source="github")
}
