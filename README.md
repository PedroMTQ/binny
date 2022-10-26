[![DOI](https://zenodo.org/badge/327396590.svg)](https://zenodo.org/badge/latestdoi/327396590)



# binny

## Requirements
At the moment binny only runs on Linux. \
Conda (and optionally, recommended Mamba) as well as Git need to be available.

## Quickstart
Here is a quick guide on the installation and test run of binny. Please check out the longer description below to set up binny on a cluster environment.

1) Clone this repository with git
```
# git clone https://github.com/a-h-b/binny.git
git clone -b sk_v17 https://github.com/ohickl/binny.git
cd binny
```

2) Create the binny environment (mamba is recommended for speed over conda).
```
# Choose env manager
my_env_manager='mamba' # or 'conda'

# Optional:
my_conda_env_path="absolute/path/to/conda/env/dir" #adjust path here

# If the conda channel priority is set to 'strict', the env creation will likely fail
# so you might need to use:
# conda config --set channel_priority flexible

${my_env_manager} env create --file workflow/envs/binny.yaml

# or
# ${my_env_manager} env create --file workflow/envs/binny.yaml --prefix ${my_conda_env_path}

# If necessary:
# conda config --set channel_priority strict
```

3) Database and Mantis setup with test run.
```
${my_env_manager} activate binny # or ${my_env_manager} activate ${my_conda_env_path}/binny

./binny --outputdir test_output \
        --assembly test/contigs_4bins.fa \
        --bam test/reads_4bins*.bam \
        --threads 4
```

If all goes well, binny will run in the current session, load the CheckM data, setup Mantis, and make and fill a directory called `test_output`. A completed run should contain four fasta files with one bin each in `test_output/bins`. 

To view all available parameters and accompanying explanations use `./binny --help`


### CheckM databases

The marker gene data file `checkm_data_2015_01_16.tar.gz` is downloaded from [here](https://data.ace.uq.edu.au/public/CheckM_databases), and the following files are processed:
* taxon_marker_sets.tsv
* tigrfam2pfam.tsv
* checkm.hmm

The processed marker gene file, `taxon_marker_sets_lineage_sorted.tsv`, can be found in the `database` directory by default and is generated using `remove_unused_checkm_hmm_profiles.py` found under
`workflow/scripts`.