## Repository notes

This repository contains all of the code required to recreate the analyses presented the manuscript, *Joint effects of host genotype and species arrival order determine plant microbiome composition and function*, ***Current Biology***, 2020. [pdf](https://dleopold.github.io/files/pubs/Leopold,%20Busby%20-%202020%20-%20Host%20Genotype%20and%20Colonist%20Arrival%20Order%20Jointly%20Govern%20Plant%20Microbiome%20Composition%20and%20Function.pdf)

Authors: Devin R. Leopold & Posy E. Busby



Assuming all dependencies are available (see below), the entire workflow can be recreated by running recipies in the makefile, which will run the scripts in the `code/` folder in the proper order, writing all output to `output/`. The processed output of the bioinformatics processing of the raw Illumina marker gene sequencing (fungal ITS) is included, so it is also possible to skip the bioinformatic processing of the raw sequencing reads and skip directly to the data analysis. 

*****

#### Data processing

  * [demux.config.json](https://github.com/dleopold/Populus_priorityEffects/blob/master/code/demux.config.json)
    * Define instructions for demultiplexing with Pheniqs.
  * [trim.sh](https://github.com/dleopold/Populus_priorityEffects/blob/master/code/trim.sh)
    * Process demultiplexed Illumina data to trim gene primers and read-through contamination.
  * [compile.R](https://github.com/dleopold/Populus_priorityEffects/blob/master/code/compile.R)
    * Process denoised Illumina data to prepare for analysis. Includes removing host contamination, collapsing denoised sequence variants to 99% OTUs, identifying focal taxa, and removing samples with poor coverage.
  * [denoise.R](https://github.com/dleopold/Populus_priorityEffects/blob/master/code/denoise.R)
    * Denoise Illumina data to identify amplicon sequence variants with DADA2.

#### Analysis

  * [biasEstimates.R](https://github.com/dleopold/Populus_priorityEffects/blob/master/code/biasEstimates.R)
    * Use mock community data to estimate taxon specific biases in the Illumina sequence data.
  * [jsdModels.R](https://github.com/dleopold/Populus_priorityEffects/blob/master/code/jsdModels.R)
    * Multvariate test using joint species distribution models with mvabund.
  * [priorityEffects.R](https://github.com/dleopold/Populus_priorityEffects/blob/master/code/priorityEffects.R)
    * Explore the benefit of preemptive colonization for the fungal species used as early colonists in the immigration history treatments.
  * [rustAnalyses.R](https://github.com/dleopold/Populus_priorityEffects/blob/master/code/rustAnalyses.R)
    * Test effects of experimental treatments on leaf rust severity and make corresponding figures.
  * [rustCor.R](https://github.com/dleopold/Populus_priorityEffects/blob/master/code/rustCor.R)
    * Explore possible correlations between relative abundance of foliar fungi and rust severity.
  * [rustSusceptibility.R](https://github.com/dleopold/Populus_priorityEffects/blob/master/code/rustSusceptibility.R)
    * Look at baseline rust susceptibility in uninoculated plants.
    
#### Figures
  * [communityFigure.R](https://github.com/dleopold/Populus_priorityEffects/blob/master/code/communityFigure.R)
    * Make the multi-panel figure of variation in community composition.
  * [mapFigS1.R](https://github.com/dleopold/Populus_priorityEffects/blob/master/code/mapFigS1.R) 
    * make map showing geographic origins of *P. trichocarpa* genotypes.

#### Other 

  * [colors.R](https://github.com/dleopold/Populus_priorityEffects/blob/master/code/colors.R)
    * Define color palettes used in figures.
  * [Rfunctions.R](https://github.com/dleopold/Populus_priorityEffects/blob/master/code/Rfunctions.R)
    * Some custom R functions used by other scripts.

*****

#### Dependencies

Larger data files are not included in the repository and need to be downloaded in order to reproduce the bioinformatics workflow. 

* Raw MiSeq data
  * Download the raw MiSeq data in fastq format from the NCBI Sequence Read Archive, BioProject ID [PRJNA605581](https://www.ncbi.nlm.nih.gov/bioproject/605581)
  * The fastq files should be placed in a folder in the project directory, `output/demux/`.
Because the archived files are demultiplexed, the initial demultiplexing step in the makefile, `make demux`, is not necessary.

* Taxonomic assignment of non-target taxa
  * Download the UNITE v8.2 databast [here](https://dx.doi.org/10.15156/BIO/786372).
  * Move to `data/referenceDBs/sh_general_release_dynamic_s_04.02.2020.fasta.gz`.
  
* Populus trichocarpa genome
  * Download the v3 genome [here](https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_Ptrichocarpa_er).
  * Create an searchable database, using `bowtie2-build`, and move the indexed files to `data/referenceDBs/Ptri.v.3.db`.

* Software dependencies (versions used).
  * bowtie2 (v2.3.5)
  * cutadapt (v2.7)
  * pheniqs (v2.0.4)
  * R (3.6.2)
  * SeqPurge (v2019_11)

* R packages
  * betareg (v3.1-3)
  * Biostrings (v2.54.0)
  * cowplot (v1.0.0)
  * dada2 (v1.15)
  * DECIPHER (v2.14)
  * doMC (v1.3.6)
  * emmeans (v1.4.4)
  * foreach (v1.4.8)
  * ggbeeswarm (v0.6.0)
  * ggtext (v0.1.0)
  * ggnewscale (v0.4.1)
  * ggthemes (v4.2.0)
  * ggvegan (v0.1-0)
  * lmtest (v0.9-37)
  * magrittr (v1.5)
  * MASS (v7.3-51.5)
  * mvabund (v4.0.1)
  * patchwork (v1.0.0)
  * philentropy (v0.4.0)
  * phyloseq (v1.30.0)
  * scico (v1.1.0)
  * ShortRead (v1.44.3)
  * tidyverse (v1.3)
  * vegan (v2.5-6)
