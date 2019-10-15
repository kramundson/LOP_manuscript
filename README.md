# LOP_manuscript

Analyses are divided among subdirectories in ```experiments```.
Each is a Snakemake workflow that is written to be run in succession. Datasets generated
in earlier experiments are used in later experiments, so executing each in order is crucial.

## Getting started

Install miniconda3

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# follow prompts
# yes to putting conda in your path

source $HOME/.bashrc

# change location of envs_dirs
# default is home directory but conda crashes if home folder is not writable
conda config --add envs_dirs ./.conda/envs

# change location of pkgs_dirs
# default is home directory but conda crashes if home folder is not writable
conda config --add pkgs_dirs ./.conda/pkgs
```

Build conda environment

```
conda env create -f environment.yaml
# follow prompts to finish environment build

# activate environment
conda activate LOP_manuscript
```

Install R package MeanShift from CRAN archive

```
wget https://cran.r-project.org/src/contrib/Archive/MeanShift/MeanShift_1.1-1.tar.gz
tar -xzvf MeanShift_1.1-1.tar.gz
# todo figure out the rest of the command line install on server
```

## Experiments

0. Build reference genome

```
cd experiments/ref
snakemake --cores 8 > ref-build.out 2> ref-build.err
```

1. Read processing, alignment, and joint variant calling of Alca Tarma, haploid inducers
PL4 and IVP101, and three dihaploids of Alca Tarma

```
cd ../1_high_cov_variants
snakemake -s 1_init_files.snakes > 1_smk.out 2> 1_smk.err
snakemake -s 2_fastq_vcf.snakes --cores 64 > 2_smk.out 2> 2_smk.err
```
 
2. Dosage analysis of 167 dihaploids sequenced at low coverage

```
cd ../2_low_cov_cnv
snakemake --cores 64 > smk.out 2> smk.err
```
 
3. Generate fake hybrids for power test of dosage variants in LOP dihaploid population

```
cd ../3_fake_hybrids
snakemake -s 1_make_simhybs.snakes --cores 64 > 1_smk.out 2> 1_smk.err
snakemake -s 2_fastq_to_bam.snakes --cores 64 > 2_smk.out 2> 2_smk.err
```

4. Binned SNP genotype analysis of 167 dihaploids. Uses alignments from Experiment 2.

```
cd ../4_low_cov_snp
snakemake --cores 24 > smk.out 2> smk.err
```

## Results

Figures and datasets described in the manuscript, as well as datasets used for analyses
generated by the experiments above are deposited in the ```results``` folder.

For experiment details, see the ```README.md``` file located in each experiment subdirectory.