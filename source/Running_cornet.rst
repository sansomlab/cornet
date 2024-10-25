Running Cornet 
==================================

You are recommended to read the `WGCNA FAQ <https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?rlkey=thqg8wlpdn4spu3ihjuc1kmlu&e=1&dl=0/>`_ before starting analysis.

The following command loads the pipeline.yml configuration file:

.. code-block:: Bash

    python ~/path/to/cornet/pipelines/pipeline_wgcna.py config

The key parameters for the analysis are defined in this pipeline.yml file - these will need to be set on a case by case basis. The input data parameters in pipeline.yml should be set as descibed on the :ref:`Data Preparation` page.

In general use of a signed or signed hybrid network and bicor correlation is recommended (see cautionary notes and advice in the `WGCNA FAQ <https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?rlkey=thqg8wlpdn4spu3ihjuc1kmlu&e=1&dl=0/>`_).

Tasks can be run by using an appropriate task name in the command:

 .. code-block:: Bash
    
    python ~/path/to/cornet/pipelines/pipeline_wgcna.py make __task_name__

The cleanData, softPower and detectModules tasks should be run to check parameterisation before running the full pipeline. The full pipeline can then be run using the 'make full' command. 


.. note:: 
    If you need to rerun any tasks of the pipeline then the appropriate sentinel file needs to be deleted prior to rerunning the task. 

Data Cleaning 
-------------
After running "make cleanData" inspect the plots in "wgcna.dir/clean.dir" to check the thresholds for sample exclusion etc are appropriate (see the "clean" section of pipeline.yml). If you need to update parameters "rm" the "wgcna.dir/modules.dir" (or delete the sentinel file) before re-running the task.

The parameters in the clean section of pipeline.yml refer to arguments for the `goodSampleGenes <https://www.rdocumentation.org/packages/WGCNA/versions/1.72-5/topics/goodSamplesGenes>`_ function from the R WGCNA package. The cut_height parameter is the cut height used to exclude samples from analysis, as shown in the 'sampleClustering' dendrogram output. The 'sampleClustering_with_traits' output can be inspected to assess the clustering of samples by trait. 

Soft Power selection
--------------------
After running "make softPower" it is essential to inspect the network toplogy plots in the "wgcna.dir/soft.power.dir" and to update the yml with an appropriate soft power as guided by the WGCNA documentation and tutorials before continuing (See `WGCNA tutorial <https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?rlkey=thqg8wlpdn4spu3ihjuc1kmlu&e=2&dl=0%2F.>`_). 

Module detection
----------------
After running "make detectModules" inspect the plots in the "wgcna.dir/modules.dir" to check that you are happy with the module merging and update the pipeline.yml "module" section accordingly. If you need to update parameters "rm" the "wgcna.dir/modules.dir" (or delete the sentinel file) before re-running the task.

This task should be run to check module merging before running the full pipeline. In particular the min_size parameter dictates the minimumn number of genes included in a module, the deepsplit parameter provides control over the cluster separation sensitivity, and the diss_threshold dictates the dissimilarity threshold used to merge modules, this is visualised on the 'Clustering of module eigengenes' plot output. 

In general we use a 'signed hybrid' network and 'bicor' correlation unless any of the assumptions described in the `WGCNA FAQ <https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?rlkey=thqg8wlpdn4spu3ihjuc1kmlu&e=1&dl=0/>`_ are violated.

.. note:: 
    The blockwise module detection functionality is experimental and untested.

Full run
--------
Once you are happy with the parameterisation of the previous steps you can build a folder containing a PDF report, module membership tsv file, geneset xlsx document, eigengene matrix tsv and eigengene expression heatmap by calling the pipeline with "make report".

The pipeline performs testing for gene set over representation in modules using `gsfisher <https://github.com/sansomlab/gsfisher>`_, generating both summary plots and a full table of results.

Additional note
----------------
Support for blockwise detection is experimental and untested. In general we have not needed yet to use blockwise detection (but we are fortunate to have access to well-resourced cluster nodes).