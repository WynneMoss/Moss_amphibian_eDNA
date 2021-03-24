# Moss et al. (in review)
Repository contains scripts and raw data for a study of the relative efficacy of environmental DNA techniques for amphibian monitoring.

Contents are for *in prep* manuscript: Moss *et al.* Navigating the tradeoffs between environmental DNA and conventional field surveys for improved amphibian monitoring.

Permanently archived at: [Pending]

##  Contents:
(1) Curated reference databases used in metabarcoding analysis (GenBank/fasta format) [(here)](https://github.com/WynneMoss/Moss_amphibian_eDNA/tree/main/01_reference_database)

(2) Notebook to run metaBEAT pipeline for processing Illumina metabarcoding data
[(here)](https://github.com/WynneMoss/Moss_amphibian_eDNA/tree/main/02_metaBEAT)

(3) NCBI Sequence read archive (SRA) accession numbers for raw Illumina metabarcoding data 
[(here)](https://github.com/WynneMoss/Moss_amphibian_eDNA/tree/main/03_raw_reads)

(4) Taxonomic assignment results (metabarcoding)  [(here)](https://github.com/WynneMoss/Moss_amphibian_eDNA/tree/main/04_taxonomic_assignment)

(5) Data - both raw and reformatted. 
[(here)](https://github.com/WynneMoss/Moss_amphibian_eDNA/tree/main/05_data)

Includes amphibian detection data from field, qPCR, and metabarcoding as well as survey covariates.

(6) Scripts for data re-formatting and analysis in R
[(here)](https://github.com/WynneMoss/Moss_amphibian_eDNA/tree/main/06_scripts)

Includes scripts for modifying and reformatting existing code into aggregated datasets and occupancy model objects.
Also includes scripts for running statistical analysis and all plots for manuscript.


## Instructions to set up dependencies for data processing and analyses

To facilitate full reproducibility of our analyses, we provide Jupyter notebooks illustrating our workflow and all necessary associated data in this repository.

Illumina data was processed (from raw reads to taxonomic assignment) using the [metaBEAT](https://github.com/HullUni-bioinformatics/metaBEAT) pipeline. The pipeline relies on a range of open bioinformatics tools, which we have wrapped up in a self-contained docker image that includes all necessary dependencies [here](https://hub.docker.com/r/chrishah/metabeat/).


## Setting up the environment

In order to retrieve scripts and associated data (reference sequences, sample metadata etc.), start by cloning this repository to your current directory:

```
git clone --recursive https://github.com/WynneMoss/Moss_amphibian_eDNA.git
```

In order to make use of our self contained analysis environment, you will have to install Docker on your computer. Docker is compatible with all major operating systems, but see the Docker documentation for details. On Ubuntu, installing Docker should be as easy as:

```
sudo apt-get install docker.io
```

Once Docker is installed, you can enter the environment by typing:

```
sudo docker run -i -t --net=host --name metaBEAT -v $(pwd):/home/working chrishah/metabeat /bin/bash
```

This will download the metaBEAT image (if not yet present on your computer) and enter the 'container' i.e. the self contained environment (**NB:** ```sudo``` may be necessary in some cases). With the above command, the container's directory ```/home/working``` will be mounted to your current working directory (as instructed by ```$(pwd)```). In other words, anything you do in the container's ```/home/working``` directory will be synced with your current working directory on your local machine.


## Data processing workflow as Jupyter notebooks

Raw illumina data will be deposited on the NCBI SRA:
- Study: 
- BioProject: 
- BioSample accessions: 
- SRA accessions: 


The sample specific accessions can be found [here](https://github.com/WynneMoss/Moss_amphibian_eDNA/tree/main/03_raw_reads/Sample_accessions.tsv). Before following the workflow for data processing, you'll need to download the raw reads from the SRA. To download the raw read data, you can follow the steps in this [Jupyter notebook](https://github.com/WynneMoss/Moss_amphibian_eDNA/tree/main/tree/main/03_raw_reads/How_to_download_from_SRA.ipynb).

With the data in place, you should be able to fully reproduce our analyses by following the steps outlined in the [Jupyter notebook](https://github.com/WynneMoss/Moss_amphibian_eDNA/tree/main/tree/main/02_metaBeat/CA_2018_pond_eDNA_metabarcoding_analysis.ipynb).

The workflow illustrated in the notebooks assumes that the raw Illumina data is present in a directory ```raw_reads``` at the base of the repository structure and that the files are named according to the following convention: 'sampleID-marker', followed by '_R1' or '_R2' to identify the forward/reverse read file respectively. SampleID must correspond to the first column in the file ```Sample_accessions.tsv``` [here](https://github.com/WynneMoss/Moss_amphibian_eDNA/tree/main/03_raw_reads/Sample_accessions.tsv).


