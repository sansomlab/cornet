Running Cornet 
==================================

You are recommended to read the `WGCNA FAQ <https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?rlkey=thqg8wlpdn4spu3ihjuc1kmlu&e=1&dl=0/>`_ before starting analysis.

The pipeline.yml file can be accessed by running the following command:

.. code-block:: Bash

    python ~/path/to/cornet/pipelines/pipeline_wgcna.py config

The initial parameters in pipeline.yml should be set as descibed on the :ref:`Data Preparation` page.

Tasks can be run using appropriate task name in the command:

 .. code-block:: Bash
    
    python ~/path/to/cornet/pipelines/pipeline_wgcna.py make __task_name__

The cleanData, softPower and detectModules tasks should initially be run to check parameterisation before running the full pipeline. The full pipeline can be run using the 'make full' command. 


.. note:: 
    If you need to rerun any tasks of the pipeline then the appropriate sentinel file needs to be deleted prior to rerunning the task. 

Data Cleaning 
-------------
The parameters in the clean section of pipeline.yml refer to arguments for the `goodSampleGenes <https://www.rdocumentation.org/packages/WGCNA/versions/1.72-5/topics/goodSamplesGenes>`_ function from the R WGCNA package. The cut_height parameter is the cut height used to exclude samples from analysis as shown in the 'sampleClustering' dendrogram output. The 'sampleClustering_with_traits' output can be inspected assess the clustering of samples by trait. 

Soft Power selection
--------------------
This task returns network topology plots which should be inspected to inform soft power selection for downstream processing (See `WGCNA tutorial <https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?rlkey=thqg8wlpdn4spu3ihjuc1kmlu&e=2&dl=0%2F.>`_). The soft_power parameter in the pipeline.yml file should be updated to the selected soft power after running this task. 

Module detection
----------------
This task should be run to check module merging before running the full pipeline. In particular the min_size parameter dictates the minimumn number of genes included in a module, the deepsplit parameter provides control over the cluster separation sensitivity, and the diss_threshold dictates the dissimilarity threshold used to merge modules, this is visualised on the 'Clustering of module eigengenes' plot output. 

In general we use a 'signed hybrid' network and 'bicor' correlation unless any of the assumptions described in the `WGCNA FAQ <https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?rlkey=thqg8wlpdn4spu3ihjuc1kmlu&e=1&dl=0/>`_ are violated.

.. note:: 
    The blockwise module detection functionality is experimental and untested.

Full run
--------
Other outputs of the pipeline include a heatmap of module-trait relationships, a heatmap of eigengenes per sample, a plot of genes that correlate with eigenegenes, a module gene membership table (formed from genes which significantly correlate with the module eigengene) and a module eigenegene expression matrix. 

The pipeline performs testing for gene set over representation in modules using `gsfisher <https://github.com/sansomlab/gsfisher>`_ generating both summary plots and a full table of results.

The pipeline will also produce a full pdf report of the key results from this analysis, this can be found in the 'report' folder.   

