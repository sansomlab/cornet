# cornet

Pipelines for building correlation networks from gene expression data.

Currently there is one pipeline, pipeline_wgcna.py. The code has been used by the group for several projects but is still undergoing development and testing and should be considered "beta" - please do your own sanity checks if you use it!


## Installation

1. Install the cgat-core pipeline system following the instructions here: https://github.com/cgat-developers/cgat-core/

2. Clone the cornet repository, e.g.
```
   $> git clone https://github.com/sansomlab/cornet.git
```

3. In the same virtual or conda environment as cgat-core install the required python packages:
```
   $> pip install -r cornet/python/requirements.txt
```

4. To install the required R packages (R>4.0.0, bioconductor and devtools are prerequisite):
```
   $> Rscript cornet/R/install.packages.R
```

## Running the pipelines

* Documentation for pipeline_wgcna.py is [avaliable here](docs/pipeline_wgcna.md).
