## Repository notes

This repository contains all of the code required to recreate the analyses presented in a scientific manuscript, *Joint effects of host genotype and species arrival order determine plant microbiome composition and function*, currently in review.

Assuming all dependencies are available (see below), typing `make analysis` will recreate the entire workflow by running the scripts in the `code/` folder in the proper order. The processed output of the bioinformatics processing of the raw Illumina marker gene sequencing (fungal ITS) is included, so it is possible to skip the bioinformatic processing of the raw sequencing reads and skip directly to the data analysis, e.g. `make analysis`. 

### Dependencies

Larger data files are not included in the repository and need to be downloaded.

* Raw MiSeq data
  * Download the raw MiSeq data in fastq format from the NCBI Sequence Read Archive, BioProjuect ID [PRJNA605581](https://www.ncbi.nlm.nih.gov/bioproject/605581)
  * The fastq files should be placed in a folder in the project directory, `output/demux/`.
Because the archived files are demultiplexed, the initial demultiplexing step in the makefile, `make demux`, is not necessary.

* Taxonomic assignment of non-target taxa
  * Download the UNITE v8.2 databast [here](https://dx.doi.org/10.15156/BIO/786372).
  * Move to `data/referenceDBs/sh_general_release_dynamic_s_04.02.2020.fasta.gz`.
  
* Populus trichocarpa genome
  * Download the v3 genome [here](https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_Ptrichocarpa_er).
  * Create an searchable database, using `bowtie2-build`, and move the indexed files to `data/referenceDBs/Ptri.v.3.db`.

* Software dependencies can be found in the individual scripts in the `code/` folder.





