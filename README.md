# Lockpick
Analyzing Locksmith-designed probe capture methylation sequencing results
## Installing Lockpick
Clone the repository with the following command:
```
$ git clone https://github.com/RoySimons96/Lockpick.git
```
Enter the repository by `cd Lockpick`.

## Dependencies


### Set up conda environment
- Download [Conda](https://www.anaconda.com/products/individual) with Python 3.9

```
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
```

Install Conda
```
bash Anaconda3-2021.11-Linux-x86_64.sh
```

- [Snakemake](https://snakemake.readthedocs.io/) (at least v4.3.1) and [biopython-1.79](https://biopython.org/docs/1.79/api/Bio.html)

Install Mamba into your Conda-based python distribution
```
conda install -n base -c conda-forge mamba
```
Activate the Conda base environment (which now includes Mamba).
```
conda activate base
```
Create a new conda environment called ```Lockpick``` with snakemake in it.
```
mamba create -c conda-forge -c bioconda -n snakemake Lockpick
```
Activate the ```Lockpick``` conda environment.
```
conda activate Lockpick
```
Check whether Snakemake is successfully installed by running the following command:
```
snakemake --help
```
Install the following packages: biopython=1.79 .
```
conda install python=3.9
conda install -c conda-forge biopython=1.79
conda install -c bioconda cutadapt=4.5
conda install -c bioconda multiqc=1.18
conda install -c bioconda fastqc=0.12.1

```
## Running Lockpick
Go into the Lockpick directory by:
```
cd <Path to Lockpick>
```
Make sure to activate the conda environment by: 
```
conda activate Lockpick
```
```
snakemake --cores [amount of cores] --configfile [path to config file]
```

### Configuration
The configuration file which is needed to run Lockpick includes all initial parameters needed to analyze sequencing data created by capture using a Locksmith probe panel.

"path_to_config_file": The path in which the configuration file can be found.

"path_to_work_dir": The path where temporary files are created.

"path_to_scripts": The directory in which the Lockmith python scripts are situated.

"output_directory":  The path of the directory in which the output files will be created.

"chosen_panel_csv": The Locksmith designed probe panel file, used for recognizing probe arms and to construct the amplicon reference file.

"reference_genome_path": The path to the FASTA format of the human reference genome to be used.

"fastq_folder": The path to the folder in which the to be analyzed FASTQ files are located.

"target_cpg_bed": The path to the bedfile which states the targets used for probe design.

"bismark_flags": Flags used when running bismark.

"cutadapt_flags": Flags used when running cutadapt for adapter removal. 

"samples": a list of sample names to be used for analysis. 

Note that the expected format for the fastq files is {sample}_{read}.fq, where {sample} is taken from the "samples" input and {read} either denotes a 1 or a 2, for the forward and reverse read respectively.
