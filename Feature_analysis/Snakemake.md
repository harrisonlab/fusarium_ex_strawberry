# SNAKEMAKE
Make conda env and install as follows

Install mamba in base envs
    conda install -n base -c conda-forge mamba
    conda activate base

Create sankemake env
    mamba create -c conda-forge -c bioconda -n snakemake snakemake

After installing mamba, create snakmake env
    conda activate snakemake

Copy METEORE data using command below
    git clone https://github.com/comprna/METEORE.git
    cd METEORE/

Run the pipeline with your own dataset by replacing example folder in the data directory with your folder containing the fast5 files. You will use the fast5 folder name to specify your target output file in the Snakemake pipeline

    git clone https://github.com/comprna/METEORE.git
    cd METEORE/
#
