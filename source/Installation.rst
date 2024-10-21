Installation 
==================================

Dependencies
------------

Core dependencies include:

- Python3
- The cgat-core pipeline framework
- Python packages as per cornet/python/requirements.txt
- R >= 4.0
- Various R libraries (bioconductor, devtools and see cornet/R/install.packages.R)

Installation 
------------

#. Install the cgat-core pipeline system following the instructions here `https://github.com/cgat-developers/cgat-core/ <https://github.com/cgat-developers/cgat-core/>`_.

#. Clone and install the cornet repository e.g.

    .. code-block:: Bash
     
         git clone https://github.com/sansomlab/cornet.git

#. In the same virtual or conda environment as cgat-core install the required python packages::

    pip install -r cornet/python/requirements.txt

#. Make sure you have the "devtools" and "BiocManager" R libraries pre-installed by running the following command in an R shell::

    install.packages(c("devtools","BiocManager"))

#. Install the required R packages by running the following command from the bash shell::

    Rscript cornet/R/install.packages.R 