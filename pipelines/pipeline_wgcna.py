##############################################################################
#
#   Kennedy Institute of Rheumatology
#
#   $Id$
#
#   Copyright (C) 2018 Stephen Sansom
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################

"""===========================
Pipeline WGCNA
===========================

:Author: Sansom lab
:Release: $Id$
:Date: |today|
:Tags: Python

Overview
========

This pipeline wraps the WGCNA method using a set of Rscripts.


Usage
=====


Configuration
-------------

The pipeline requires a configured :file:`pipeline.yml` file.

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_seurat.py config


Input files
-----------



Dependencies
------------

This pipeline requires:

* cgat-core: https://github.com/cgat-developers/cgat-core
* R & various packages.
* Latex.


Pipeline output
===============


"""

from ruffus import *
from pathlib import Path
import sys
import os
import shutil
import glob
import sqlite3
import numpy as np
import pandas as pd
import textwrap
import subprocess
from scipy.stats.mstats import gmean
import cgatcore.experiment as E
from cgatcore import pipeline as P
import cgatcore.iotools as IOTools

# import local pipeline utility functions
from pipeline_utils import templates


# -------------------------- < parse parameters > --------------------------- #

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

# set the location of the tenx code directory
if "wgnca_dir" not in PARAMS.keys():
    PARAMS["wgcna_dir"] = Path(__file__).parents[1]
else:
    raise ValueError("Could not set the location of the tenx code directory")


# ----------------------- < pipeline configuration > ------------------------ #

# handle pipeline configuration
if len(sys.argv) > 1:
        if(sys.argv[1] == "config") and __name__ == "__main__":
                    sys.exit(P.main(sys.argv))


# ########################################################################### #
# ################# Sanity check and parse the input files ################## #
# ########################################################################### #

global EXPRESSION_DATA_PATH
EXPRESSION_DATA_PATH = PARAMS["input_expression_data"]

if not os.path.exists(EXPRESSION_DATA_PATH):
    raise ValueError("input expression data file not found")

# optional trait and metadata inputs

global TRAIT_DATA_STAT
if not PARAMS["input_trait_data"] is None:
    trait_data_path = PARAMS["input_trait_data"]

    if not os.path.isfile(trait_data_path):
        raise ValueError("trait file not found")

    else:
        TRAIT_DATA_STAT = "--traitdata=%(trait_data_path)s" % locals()

else:
    TRAIT_DATA_STAT = ""

global META_DATA_STAT
if not PARAMS["input_meta_data"] is None:
    meta_data_path = PARAMS["input_meta_data"]

    if not os.path.isfile(meta_data_path):
        raise ValueError("meta file not found")

    else:
        META_DATA_STAT = "--metadata=%(meta_data_path)s" % locals()

else:
    META_DATA_STAT = ""





@follows(mkdir("annotation.dir"))
@files(None, "annotation.dir/genesets.sentinel")
def getGenesetAnnotations(infile, outfile):
    '''Get mappings between Ensembl gene_ids and (i) Entrez ids
       and (ii) KEGG pathways.
    '''

    outdir = os.path.dirname(outfile)

    log_file = outfile.replace(".sentinel", ".log")

    statement = '''Rscript %(wgcna_dir)s/R/wgcna_fetch_geneset_annotations.R
                 --ensemblversion=%(annotation_ensembl_release)s
                 --ensemblhost=%(annotation_ensembl_host)s
                 --species=%(annotation_species)s
                 --outdir=%(outdir)s
                 &> %(log_file)s
              ''' % dict(PARAMS, **locals())

    # requires internet connectivity.
    # and the BMRC cluster is broken by to_cluster = FALSE!
    process = subprocess.Popen(statement.replace("\n", ""),
                               shell=True, stdout=subprocess.PIPE)
    process.wait()

    if process.returncode != 0:
        raise ValueError("failed to get annotation")

    IOTools.touch_file(outfile)


@files(None,
       "wgcna.dir/clean.dir/clean.sentinel")
def cleanData(infile, outfile):
    '''
    Prepare the data for a WGCNA run
    '''

    results_file = outfile.replace(".sentinel", ".RData")
    log_file = outfile.replace(".sentinel", ".log")

    #infile_path = os.path.abspath(infile)

    expression_data_path = EXPRESSION_DATA_PATH
    trait_data_stat = TRAIT_DATA_STAT

    out_dir = os.path.dirname(os.path.abspath(outfile))
    results_filename = os.path.basename(results_file)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    statement = '''Rscript %(wgcna_dir)s/R/wgcna_data_cleaning.R
                   --input=%(expression_data_path)s
                   --idcol=%(annotation_idcol)s
                   --outdir=%(out_dir)s
                   --outfilename=%(results_filename)s
                   --minfraction=%(clean_min_fraction)s
                   --minnsamples=%(clean_min_n_samples)s
                   --minngenes=%(clean_min_n_genes)s
                   --minrelativeweight=%(clean_min_relative_weight)s
                   --cutheight=%(clean_cut_height)s
                   --minsize=%(clean_min_size)s
                   %(trait_data_stat)s
                   &> %(log_file)s
                '''

    P.run(statement)

    IOTools.touch_file(outfile)


@transform(cleanData,
           regex(r"(.*)/.*/clean.sentinel"),
           r"\1/soft.power.dir/soft.power.sentinel")
def softPower(infile, outfile):
    '''
    Run the soft power analysis
    '''

    log_file = outfile.replace(".sentinel", ".log")

    clean_data = infile.replace(".sentinel", ".RData")

    job_threads = PARAMS["module_threads"]
    job_memory = PARAMS["module_memory"]

    out_dir = os.path.dirname(outfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    statement = '''Rscript %(wgcna_dir)s/R/wgcna_soft_power.R
                   --input=%(clean_data)s
                   --outdir=%(out_dir)s
                   --networktype=%(module_network_type)s
                   --adjcorfnc=%(module_adj_cor_fnc)s
                   --adjdistfnc=%(module_adj_dist_fnc)s
                   --threads=%(module_threads)s
                   &> %(log_file)s
                '''

    P.run(statement)

    IOTools.touch_file(outfile)



# ########################################################################### #
# ################### Step by Step module detection ######################### #
# ########################################################################### #


@transform(cleanData,
           regex(r"(.*)/.*/clean.sentinel"),
           r"\1/modules.dir/adjacency.sentinel")
def computeAdjacency(infile, outfile):
    '''Compute the adjacency matrix'''

    results_file = os.path.basename(outfile).replace(".sentinel", ".RData")
    log_file = outfile.replace(".sentinel", ".log")

    clean_data = infile.replace(".sentinel", ".RData")

    job_threads = PARAMS["module_threads"]
    job_memory = PARAMS["module_memory"]

    out_dir = os.path.dirname(outfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    statement = '''Rscript %(wgcna_dir)s/R/wgcna_compute_adjacency.R
                   --input=%(clean_data)s
                   --outdir=%(out_dir)s
                   --outfilename=%(results_file)s
                   --threads=%(module_threads)s
                   --softpower=%(module_soft_power)s
                   --networktype=%(module_network_type)s
                   --adjcorfnc=%(module_adj_cor_fnc)s
                   --adjdistfnc=%(module_adj_dist_fnc)s
                   &> %(log_file)s
                '''

    P.run(statement)

    IOTools.touch_file(outfile)

@transform(computeAdjacency,
           regex(r"(.*)/modules.dir/adjacency.sentinel"),
           r"\1/modules.dir/TOM.sentinel")
def computeTOM(infile, outfile):
    '''Compute the TOM'''

    results_file = os.path.basename(outfile).replace(".sentinel", ".RData")
    log_file = outfile.replace(".sentinel", ".log")

    adjacency_data = infile.replace(".sentinel", ".RData")

    job_threads = PARAMS["module_threads"]
    job_memory = PARAMS["module_memory"]

    out_dir = os.path.dirname(outfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    statement = '''Rscript %(wgcna_dir)s/R/wgcna_compute_TOM.R
                   --input=%(adjacency_data)s
                   --outdir=%(out_dir)s
                   --outfilename=%(results_file)s
                   --threads=%(module_threads)s
                   --tomtype=%(module_tom_type)s
                   &> %(log_file)s
                '''

    P.run(statement)

    IOTools.touch_file(outfile)


@transform(computeTOM,
           regex(r"(.*)/modules.dir/TOM.sentinel"),
           add_inputs(cleanData),
           r"\1/modules.dir/modules.sentinel")
def detectModules(infiles, outfile):
    '''Cluster the TOM, cut the tree and merge the modules'''

    results_file = os.path.basename(outfile).replace(".sentinel", ".RData")
    log_file = outfile.replace(".sentinel", ".log")

    tom_data, clean_data = [x.replace(".sentinel", ".RData") for x in infiles]

    job_threads = PARAMS["module_threads"]
    job_memory = PARAMS["module_memory"]

    out_dir = os.path.dirname(outfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    statement = '''Rscript %(wgcna_dir)s/R/wgcna_detect_modules.R
                   --cleandata=%(clean_data)s
                   --tomdata=%(tom_data)s
                   --outdir=%(out_dir)s
                   --outfilename=%(results_file)s
                   --threads=%(module_threads)s
                   --softpower=%(module_soft_power)s
                   --minmodulesize=%(module_min_size)s
                   --medissthreshold=%(module_diss_threshold)s
                   --adjcorfnc=%(module_adj_cor_fnc)s
                   &> %(log_file)s
                '''

    P.run(statement)

    IOTools.touch_file(outfile)


# ########################################################################### #
# ################# Blockwise module detection ############################## #
# ########################################################################### #


@transform(cleanData,
           regex(r"(.*)/.*/clean.sentinel"),
           r"\1/modules.dir/modules.sentinel")
def detectModulesBlockwise(infile, outfile):
    '''
    Run the module detection
    '''

    results_file = os.path.basename(outfile).replace(".sentinel", ".RData")
    log_file = outfile.replace(".sentinel", ".log")

    clean_data = infile.replace(".sentinel", ".RData")

    job_threads = PARAMS["module_threads"]
    job_memory = PARAMS["module_memory"]

    out_dir = os.path.dirname(outfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    statement = '''Rscript %(wgcna_dir)s/R/wgcna_detect_modules_blockwise.R
                   --input=%(clean_data)s
                   --outdir=%(out_dir)s
                   --outfilename=%(results_file)s
                   --threads=%(module_threads)s
                   --maxblocksize=%(module_block_size)s
                   --softpower=%(module_soft_power)s
                   --networktype=%(module_network_type)s
                   --adjcorfnc=%(module_adj_cor_fnc)s
                   --adjdistfnc=%(module_adj_dist_fnc)s
                   --tomtype=%(module_tom_type)s
                   --minmodulesize=%(module_min_size)s
                   --medissthreshold=%(module_diss_threshold)s
                   &> %(log_file)s
                '''

    P.run(statement)

    IOTools.touch_file(outfile)



# ########################################################################### #
# ################## integrate module detection options ##################### #
# ########################################################################### #

if PARAMS["module_detection"] == "stepwise":
    collectModules = detectModules

elif PARAMS["module_detection"] == "blockwise":
    collectModules = detectModulesBlockwise

else:
    raise ValueError('Module detection must be set to either "stepwise" or "blockwise"')


# ########################################################################### #
# #################### Module characterisation ############################## #
# ########################################################################### #


@transform(collectModules,
           regex(r"(.*)/.*/modules.sentinel"),
           add_inputs(getGenesetAnnotations, cleanData),
           r"\1/membership.dir/membership.sentinel")
def characteriseModules(infiles, outfile):
    '''
    Run the module detection
    '''

    modulesx, annotationsx, cleanx = infiles

    results_file = os.path.basename(outfile).replace(".sentinel", ".tsv")
    log_file = outfile.replace(".sentinel", ".log")

    clean_data = cleanx.replace(".sentinel", ".RData")
    module_data = modulesx.replace(".sentinel", ".RData")
    annotation_file = annotationsx.replace("genesets.sentinel", "ensembl.to.entrez.tsv.gz")

    job_threads = PARAMS["module_threads"]
    job_memory = PARAMS["module_memory"]

    out_dir = os.path.dirname(outfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    statement = '''Rscript %(wgcna_dir)s/R/wgcna_modules_vs_traits.R
                   --input=%(clean_data)s
                   --modules=%(module_data)s
                   --annotation=%(annotation_file)s
                   --idcol=%(annotation_idcol)s
                   --namecol=%(annotation_namecol)s
                   --outdir=%(out_dir)s
                   --outfilename=%(results_file)s
                   --threads=%(module_threads)s
                   &> %(log_file)s
                '''

    P.run(statement)

    IOTools.touch_file(outfile)


@transform(characteriseModules,
           regex(r"(.*)/.*/membership.sentinel"),
           r"\1/eigengenes.dir/eigengenes.sentinel")
def characteriseEigengenes(infile, outfile):
    '''
    Characterise the eigen genes
    '''

    modulesx = infile

    eigengene_path = os.path.join(os.path.dirname(modulesx),
                                  "eigengenes.tsv")

    membership_path = os.path.join(os.path.dirname(modulesx),
                                  "membership.tsv")

    log_file = outfile.replace(".sentinel", ".log")

    trait_data_stat = TRAIT_DATA_STAT
    meta_data_stat = META_DATA_STAT

    out_dir = os.path.dirname(outfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    statement = '''Rscript %(wgcna_dir)s/R/wgcna_characterise_eigengenes.R
                   --eigengenes=%(eigengene_path)s
                   --namecol=%(annotation_namecol)s
                   --membership=%(membership_path)s
                   --params=pipeline.yml
                   %(trait_data_stat)s
                   %(meta_data_stat)s
                   --figwidth=%(plot_eigengene_heatmap_width)s
                   --figheight=%(plot_eigengene_heatmap_height)s
                   --outdir=%(out_dir)s
                   &> %(log_file)s
                '''

    P.run(statement)

    IOTools.touch_file(outfile)

@active_if(PARAMS["input_genelists"]!=None)
@transform(collectModules,
           regex(r"(.*)/.*/modules.sentinel"),
           add_inputs(getGenesetAnnotations, cleanData),
           r"\1/eigengenes.dir/eigengenes.vs.genelists.sentinel")
def eigengenesVsGenelists(infiles, outfile):
    '''
    Characterise the eigen genes
    '''

    modulex, annotationx, cleanx = infiles

    log_file = outfile.replace(".sentinel", ".log")

    clean_data = cleanx.replace(".sentinel", ".RData")
    module_data = modulex.replace(".sentinel", ".RData")
    annotation_file = annotationx.replace("genesets.sentinel", "ensembl.to.entrez.tsv.gz")

    out_dir = os.path.dirname(outfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    statement = '''Rscript %(wgcna_dir)s/R/wgcna_eigengenes_vs_genelists.R
                   --input=%(clean_data)s
                   --annotation=%(annotation_file)s
                   --modules=%(module_data)s
                   --genelists=%(input_genelists)s
                   --idcol=%(annotation_idcol)s
                   --namecol=%(annotation_namecol)s
                   --outdir=%(out_dir)s
                   &> %(log_file)s
                '''

    P.run(statement)

    # write the tex snippet.
    genelists = pd.read_csv(PARAMS["input_genelists"], sep="\t")
    plot_groups = set(genelists["plot_group"])

    tex_file = os.path.join(out_dir, "genelists.tex")
    with open(tex_file, "w") as tex:

        for plot_group in plot_groups:

            pg_esc = plot_group.replace('_','\_')

            heatmap_path = os.path.join(out_dir,
                                        "genelist." + plot_group)

            heatmap_fig = {"width": "1", "height": "0.9",
                           "path": heatmap_path,
                           "caption": "Heatmap of manually curated " +\
                           pg_esc + " genes"}

            tex.write(templates.subsection % {"title": pg_esc + " genes"})
            tex.write(textwrap.dedent(
                templates.figure % heatmap_fig))
            tex.write("\n")

    IOTools.touch_file(outfile)



def parseGMTs(param_keys=["gmt_pathway_files_"]):
    '''Helper function for parsing the lists of GMT files'''

    all_files = []
    all_names = []

    for param_key in param_keys:


        gmts = [x for x in PARAMS.keys()
                if x.startswith(param_key)]

        if len(gmts) > 0:
            all_files += [PARAMS[x] for x in gmts]

            all_names += [x.replace(param_key, "")
                              for x in gmts]

    if len(all_files) == 0:
        all_files = "none"
        all_names = "none"
    else:
        all_files = ",".join(all_files)
        all_names = ",".join(all_names)

    return all_names, all_files


# ########################################################################### #
# ######################### Geneset Analysis ################################ #
# ########################################################################### #

@active_if(PARAMS["run_genesets"])
@transform(characteriseModules,
           regex(r"(.*)/membership.dir/.*.sentinel"),
           add_inputs(getGenesetAnnotations),
           r"\1/genesets.dir/geneset.analysis.sentinel")
def genesetAnalysis(infiles, outfile):
    '''
    Naive geneset over-enrichment analysis of module genes.

    Testing is performed with the gsfisher package.

    GO categories and KEGG pathways are tested by default.

    Arbitrary sets of genes cat be supplied as GMT files
    (e.g. such as those from MSigDB).
    '''

    membershipx, genesetAnno = infiles

    membership_file = membershipx.replace(".sentinel", ".tsv")

    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    anno = os.path.join(os.path.dirname(genesetAnno),
                        "ensembl.to.entrez.tsv.gz")

    kegg_pathways = os.path.join(os.path.dirname(genesetAnno),
                                 "kegg_pathways.rds")

    param_keys = ["gmt_celltype_files_",
                  "gmt_pathway_files_"]
    gmt_names, gmt_files = parseGMTs(param_keys=param_keys)

    # get the set of modules
    xx = pd.read_csv(membership_file, sep="\t")
    modules = [x for x in set(xx["module"].values)]

    job_memory = "20G"

    statements = []

    species = PARAMS["annotation_species"]
    wgcna_dir = PARAMS["wgcna_dir"]
    idcol = PARAMS["annotation_idcol"]

    for module in modules:

        logfile = os.path.join(outdir,
                               "geneset.analysis." + module + ".log")

        statements.append('''Rscript %(wgcna_dir)s/R/wgcna_modules_vs_genesets.R
                            --input=%(membership_file)s
                            --module=%(module)s
                            --species=%(species)s
                            --annotation=%(anno)s
                            --idcol=%(idcol)s
                            --kegg_pathways=%(kegg_pathways)s
                            --gmt_names=%(gmt_names)s
                            --gmt_files=%(gmt_files)s
                            --outdir=%(outdir)s
                            &> %(logfile)s
                      ''' % locals())

    P.run(statements)

    IOTools.touch_file(outfile)

@active_if(PARAMS["run_genesets"])
@transform(genesetAnalysis,
           regex(r"(.*)/.*.sentinel"),
           add_inputs(characteriseModules),
           r"\1/summarise.geneset.analysis.sentinel")
def summariseGenesetAnalysis(infiles, outfile):
    '''
    Summarise the geneset over-enrichment analyses of cluster marker genes.

    Enriched pathways are summarised in an Excel table and a heatmap.
    '''

    outdir = os.path.dirname(outfile)

    infile, membershipx = infiles

    membership_file = membershipx.replace(".sentinel", ".tsv")

    # need to sort out the dependencies properly!
    genesetdir = os.path.dirname(infile)

    param_keys = ["gmt_celltype_files_",
                  "gmt_pathway_files_"]
    gmt_names, gmt_files = parseGMTs(param_keys=param_keys)

    # get the set of modules
    xx = pd.read_csv(membership_file, sep="\t")
    modules = [x for x in set(xx["module"].values)]
    module_list = ",".join(modules)

    job_memory = "20G"

    logfile = outfile.replace(".sentinel", ".log")

    use_adjusted = str(PARAMS["genesets_use_adjusted_pvalues"]).upper()
    show_common = str(PARAMS["genesets_show_common"]).upper()

    show_detailed = str(PARAMS["genesets_show_detailed"])

    statement = '''Rscript %(wgcna_dir)s/R/wgcna_summariseGenesets.R
                         --genesetdir=%(genesetdir)s
                         --gmt_names=%(gmt_names)s
                         --show_detailed=%(show_detailed)s
                         --modulelist=%(module_list)s
                         --mingenes=%(genesets_min_fg_genes)s
                         --pvaluethreshold=%(genesets_pvalue_threshold)s
                         --padjustmethod=%(genesets_padjust_method)s
                         --useadjusted=%(use_adjusted)s
                         --minoddsratio=%(genesets_min_odds_ratio)s
                         --showcommon=%(show_common)s
                         --outprefix=%(outdir)s/cluster.genesets
                         --prefix=genesets
                         --plotdirvar=clusterGenesetsDir
                    &> %(logfile)s
                      '''
    P.run(statement)

    IOTools.touch_file(outfile)


# ------------------- < within cluster geneset analysis > ------------------- #

@transform(summariseGenesetAnalysis,
           regex("(.*)/genesets.dir/summarise.geneset.analysis.sentinel"),
           add_inputs(softPower, characteriseModules,
                      characteriseEigengenes, eigengenesVsGenelists),
           r"\1/latex.dir/report.vars.sty")

def latexVars(infiles, outfile):
    '''
    Prepare a file containing the latex variable definitions.
    '''

    infile = infiles[0]

    rundir = Path(outfile).parents[1]

    if not os.path.exists(os.path.dirname(outfile)):
        os.mkdir(os.path.dirname(outfile))

    vars = {"wgcnaDir": PARAMS["wgcna_dir"],
            "reportTitle": PARAMS["report_title"],
            "reportAuthor": PARAMS["report_author"],
            "minFraction": PARAMS["clean_min_fraction"],
            "minSamples": PARAMS["clean_min_n_samples"],
            "minGenes": PARAMS["clean_min_n_genes"],
            "cutHeight": PARAMS["clean_cut_height"],
            "clusterMinSize": PARAMS["clean_min_size"],
            "softPower": PARAMS["module_soft_power"],
            "detection": PARAMS["module_detection"],
            "networkType": PARAMS["module_network_type"].replace("_","-"),
            "adjCorFunction": PARAMS["module_adj_cor_fnc"],
            "adjDistFunction": PARAMS["module_adj_dist_fnc"],
            "tomType": PARAMS["module_tom_type"],
            "minSize": PARAMS["module_min_size"],
            "dissThreshold": PARAMS["module_diss_threshold"],
            "cleanDir": os.path.join(rundir,"clean.dir"),
            "powerDir": os.path.join(rundir,"soft.power.dir"),
            "membershipDir": os.path.join(rundir, "membership.dir"),
            "moduleDir": os.path.join(rundir, "modules.dir"),
            "eigengeneDir": os.path.join(rundir, "eigengenes.dir"),
            "genesetDir": os.path.join(rundir, "genesets.dir"),
            "clusterGenesetsDir": os.path.join(rundir, "genesets.dir")}

    with open(outfile, "w") as ofh:
        for command, value in vars.items():

            ofh.write("\\newcommand{\\" + command + "}{" + str(value) + "}\n")



@transform(latexVars,
           regex("(.*)/report.vars.sty"),
           r"\1/summaryReport.pdf")
def summaryReport(infile, outfile):
    '''
    Prepare a PDF summary report.
    '''

    outfile_name = os.path.basename(outfile)
    jobName = outfile_name[:-len(".pdf")]

    outdir = os.path.dirname(outfile)
    rundir = Path(outdir).parents[0]

    compilation_dir = os.path.join(outdir, ".latex_compilation.dir")

    latexVars = os.path.join(outdir, "report.vars.sty")

    try:
        shutil.rmtree(compilation_dir)
    except FileNotFoundError:
        pass

    os.mkdir(compilation_dir)

    # get the latex variables
    statement = '''pdflatex -output-directory=%(compilation_dir)s
                            -jobname=%(jobName)s
                            %(draft_mode)s
      '\\input %(latexVars)s
       \\def\\reportTitle{pipeline\\_wgcna.py: summary report}
                '''
    # get the intro
    statement += '''
      \\input %(wgcna_dir)s/pipelines/pipeline_wgcna/introReport.tex
      '''

    # add the section to visualise the soft power
    statement += '''
         \\input %(wgcna_dir)s/pipelines/pipeline_wgcna/paramSection.tex
        '''

    statement += '''
         \\input %(wgcna_dir)s/pipelines/pipeline_wgcna/cleanSection.tex
        '''

    statement += '''
         \\input %(wgcna_dir)s/pipelines/pipeline_wgcna/moduleSection.tex
        '''

    statement += '''
         \\input %(wgcna_dir)s/pipelines/pipeline_wgcna/eigengeneSection.tex
        '''

    if not PARAMS["input_genelists"] == None:
        statement += '''
        \\input %(wgcna_dir)s/pipelines/pipeline_wgcna/genelistSection.tex
        '''

    statement += '''
         \\input %(wgcna_dir)s/pipelines/pipeline_wgcna/membershipSection.tex
        '''

    if(PARAMS["run_genesets"]):
        statement += '''
        \\input %(wgcna_dir)s/pipelines/pipeline_wgcna/genesetSection.tex
                     '''

    statement += '''\\input %(wgcna_dir)s/latex/endmatter.tex'
    '''

    # Deliberately run twice - necessary for LaTeX compilation..
    draft_mode = "-draftmode"
    P.run(statement)

    draft_mode = ""
    P.run(statement)

    # Move the compiled pdfs to report.dir
    shutil.move(os.path.join(compilation_dir, outfile_name),
                outfile)


@follows(mkdir("report.dir"))
@transform(summaryReport,
           regex(r"wgcna.dir/latex.dir/summaryReport.pdf"),
           r"report.dir/report.sentinel")
def report(infile, outfile):
    '''
    Link output files to a directory in the "reports.dir" folder.

    '''

    out_dir = "report.dir"

    try:
        os.stat(out_dir)
    except FileNotFoundError:
        os.mkdir(out_dir)

    run_dir = "wgcna.dir"

    targets = {os.path.join(run_dir,"eigengenes.dir","eigengen_heatmap.png"): "module.eigengene.heatmap.png",
               os.path.join(run_dir,"genesets.dir","cluster.genesets.xlsx"): "module.genesets.xlsx",
               os.path.join(run_dir,"latex.dir","summaryReport.pdf"): "summary.report.pdf",
               os.path.join(run_dir, "membership.dir", "eigengenes.tsv"): "module.eigengene.expression.matrix.tsv",
               os.path.join(run_dir, "membership.dir", "membership.tsv"): "module.gene.membership.tsv"
               }

    for source_path, target_name in targets.items():

        if os.path.exists(source_path):

            target_path = os.path.join(out_dir, target_name)

            os.symlink(os.path.relpath(source_path, start=out_dir),
                       target_path)

    IOTools.touch_file(outfile)




# ########################################################################### #
# ##################### full target: to run all tasks ####################### #
# ########################################################################### #

@follows(report)
def full():
    pass


# ------------------- < ***** end of pipeline **** > ------------------------ #

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
