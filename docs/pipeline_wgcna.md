# Running pipeline_wgcna.py

This pipeline wraps the [WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/) R package. The original publication for WGCNA can be found [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559).

It add a few extra plots, testing for geneset over represenation in the modules with [gsfisher](https://github.com/sansomlab/gsfisher) and investigation of correlation of set of arbitrary lists of genes of interest with module eigen genes.

## General notes

* It is strongly recommended to start by reading the [step-by-step WGCNA R tutorial](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html) and the [WGCNA FAQ](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html)

* The key parameters for the analysis are exposed in the pipeline.yml file - these will need to be set on a case by case basis.

* In general use of a signed or signed hybrid network is [recommended](https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/) (we typically use a signed hybrid network as is recommended in the [WGCNA FAQ](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html)).

* In general, if the cautionary notes do not apply, we prefer to use bicor correlation as recommended in the [WGCNA FAQ](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html).


## Recommended way to run the pipeline to investigate key parameter choices

For documentation of the input files please see the "input" section of the yml file.

* As above, we tend to use a "signed hybrid network" and "bicor" correlation for gene expression analysis unless any of the assumptions (see e.g. the FAQ) are violated.

* The pipeline should be run in stages to check the data cleaning, soft power and module detection parameterisation before proceeding with a full run.

* After running "make cleanData" inspect the plots in "wgcna.dir/clean.dir" to check the threshold for sample exclusion etc are appropriate (see the "clean" section of pipeline.yml). If you need to update parameters "rm" the "wgcna.dir/modules.dir" (or delete the sentinel file) before re-running the task.

* After running "make softPower" it is essential to inspect the network toplogy plots in the "wgcna.dir/soft.power.dir" and to update the yml with an appropriate soft power as guided by the WGCNA documentation and tutorials before continuing.

* After running "make detectModules" inspect the plots in the "wgcna.dir/modules.dir" and check that you are happy with the module merging and update the pipeline.yml "module" section accordingly. If you need to update parameters "rm" the "wgcna.dir/modules.dir" (or delete the sentinel file) before re-running the task.

* Once you are happy with the parameterisation of these steps you can build a folder containing the PDF report, module membership tsv file, geneset xlsx document, eigengene matrix tsv and eigengene expression heatmap by calling the pipeline with "make report".


## Additional notes on running the pipeline

* The geneset anaysis can be toggled on or off by setting the run_genesets parameter to True|False - this is useful for initial runs where parameter choices are being explored and optimised

* Export for blockwise detection is experimental and untested. In general we have not needed to use blockwise detection (we have access to well-resourced cluster nodes).