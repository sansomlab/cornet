Data Preparation 
==================================

Expression Data
----------------

Please see the `WGCNA FAQ <https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?rlkey=thqg8wlpdn4spu3ihjuc1kmlu&e=1&dl=0/>`_ for advice on filtering out 
low expressed genes, transforming and normalising the data. 
The expression data needs to be saved as a tab-separated file. The gene idenfiers need to be in a column named with a value which is passed to the "idcol" annotation in pipeline.yml.

Trait Data
-----------
The trait data should be a tab-separated file containing quantitative trait data with the samples in a column titled 'sample_name'.

.. note:: 
    Each level of the trait data should be encoded as 0 and 1.
    
    e.g. sex should be encoded as two columns: sex_male and sex_female. Each column should contain 0 for when this is false (e.g. sex_male=0 is a female) and 1 for when this is true 
    (e.g. sex_male=1 is a male) 

Metadata
---------
The metadata should be a tab-separated file containing quantitative and/or categorical variables. The samples should be in a column titled 'sample_name'.

Gene List
----------
A tab-separated file containing a list of genes to compre to the module eigenegenes should be provided.  

This table should contain gene identifiers in a column named the same as the "idcol" annotation and the expression data, along with a "gene_name" column containing equivalent gene names, and 2 columns named "category" and "plot_group", each containing character strings, to dictate plotting behaviour.  