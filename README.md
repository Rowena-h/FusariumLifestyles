# *Fusarium* Lifestyles
 
Bioinformatics analysis pipeline for Hill et al. (in prep) Lifestyle transitions in fusarioid fungi are frequent and lack clear genomic signatures.

![Pipeline workflow](pipeline.png)

The pipeline was written for and run on Queen Mary University of London's [Apocrita HPC facility](http://doi.org/10.5281/zenodo.438045) which uses the Univa Grid Engine batch-queue system.

## 1 Assembly
`cd assembly`
### 1 Read quality control
`cd assembly/reads`
1. `qsub trimmomatic.sh` - trims raw reads using Trimmomatic.
2. `qsub fastqc.sh` - after trimming, checks read quality with FastQC.

### 2 *De novo* genome assembly
`cd assembly/denovo_assembly`
1. `./submit_assembly` - makes new directory and submits job scripts for each assembly tool (ABySS, MEGAHIT, SPAdes).
2. `./abyss_comp.sh` - after ABySS has finished running for all kmer sizes, compare the assembly stats to choose 'best' kmer size.

### 3 Polishing
`cd assembly/polishing`
`qsub polish.sh` - for each strain and assembly tool, maps raw reads to assembly and calculates mapping statistics with BWA-MEM and SAMtools and polishes the assembly with Pilon.

After completing [1 Assembly 4 Assessment](https://github.com/Rowena-h/FusariumLifestyles/tree/main/assembly/assessment) and uploading to NCBI:

`./ncbi_filter.sh` - removes sequences identified as mitochondrial or duplicates by NCBI.

### 4 Assessment
`cd assembly/assessment`
1. `./submit_assessment` - submits QUAST, BUSCO and BLAST jobs.
2. `qsub blobtools.sh` - after the BLAST has finished, submit BlobTools.

## 2 Annotation
`cd annotation`
### 1 Repeatmasking
`cd annotation/repeat_masking`
1. `qsub repeatmodeler.sh` - makes custom repeat library for each strain.
2. `qsub repeatmasker.sh` - after repeat modelling, uses custom repeat library to softmask the assembly.

### 2 MAKER pipeline
`cd annotation/maker`
#### Round 1
`cd annotation/maker/round1`
`qsub maker.sh` - first run of MAKER using ESTs and proteins (as indicated in `.ctl` files)

#### Round 2
`cd annotation/maker/round2`
1. `qsub training_snap/snap.sh` - train SNAP using gene models from the first MAKER round.
2. `qsub maker2.sh` - second run of MAKER using trained SNAP (as indicated in `.ctl` files).

#### Round 3
`cd annotation/maker/round3`
1. `qsub training_snap2/snap2.sh` - trains SNAP again using gene models from the second MAKER round.
2. `qsub maker3.sh` - third run of MAKER using second trained SNAP (as indicated in `.ctl` files).
3. `qsub rename.sh` - after obtaining unique locus tags from e.g. NCBI, renames IDs in gff and fasta files.
4. `qsub gag.sh` - removes introns <10bp, removes terminal Ns and corrects start and stop codons in gff file for NCBI compliance.

## 3 Orthology inference
`cd orthology_inference`
1. `./submit_protein_download.sh` - submits download of predicted protein sets of *Fusarium* strains from NCBI.
2. `qsub orthofinder.sh` - submits orthology inference.

## 4 Phylogenomics
`cd phylogenomics`
1. `./submit_alignment.sh` - submits alignment of single copy orthogroups from OrthoFinder with MAFFT followed by trimming with BMGE.
2. `./concat.sh` - concatenate single copy orthogroup alignments and prepare partition files.
3. `./submit_modeltestng.sh` - runs ModelTest-NG on all single copy orthogroups (in computationally tractable chunks).
4. `./submit_speciestrees_concatenation` - submits concatenation-based species tree methods (RAxML-NG and IQ-TREE).
5. `./submit_RAxML-NG_genetrees.sh` - submits RAxML-NG for individual gene trees.
6. `qsub species_tree/astral.sh` - submits coalescent-based species tree method (ASTRAL-III) using genes trees.

## 5 Divergence time estimation
`cd divergence_time_estimation`
### 1 SortaDate
`cd divergence_time_estimation/sortadata`
`qsub sortadate` - reroot gene and species trees and run with SortaDate to filter for top ten 'clock-like' genes.

### 2 MCMCTree
`cd divergence_time_estimation/mcmctree`
1. `qsub submit_mcmctree_dating_step1.sh` - submits first step of approximate likelihood divergence time estimation with protein data (see [tutorial](http://abacus.gene.ucl.ac.uk/software/MCMCtree.Tutorials.pdf).
2. `Rscript estimate_rate.r` - estimates the scaling parameter for the substitution rate prior using the species tree.
3. `./submit_mcmctree_dating_step2.sh` - submitssecond step of approximate likelihood estimation for both independent and correlated rates relaxed clock models.

## 6 CSEP prediction
`cd CSEP_prediction`
1. `./submit_CSEPprediction.sh` - submits all programmes in the CSEP prediction pipeline.
2. `./submit_CSEPfilter` - runs all programme results through CSEPfilter to produce lists of CSEPs.
3. `qsub submit_orthogroupparsing.sh` - makes presence absence matrices of orthogroups and CSEPs for all taxa and calculate stats for orthogroups.


## 7 Lifestyle comparison
`cd lifestyle_comparison`
`./submit_lifestyletest.sh` - submits lifestyle test on orthogroup and CSEP presence absence matrices; `run_edited.py` is modified from the original script found [here](https://github.com/fantin-mesny/Effect-Of-Biological-Categories-On-Genomes-Composition).

## 8 Selection
`cd selection`
1. `qsub gbff_files/ncbi_gbff_download.sh` - downloads GBFF files from NCBI; also need Ilysp transcripts in `gbff_files` directory.
2. `./submit_pal2nal.sh` - submits script to pull corresponding nucleotides for all proteins and prepares codon alignments.
3. `./submit_hyphy.sh` - submits HyPhy dN/dS selection methods.

### Codon optimisation
`cd selection/codon_optimisation`
1. `./pull_ribosomes.sh` - extracts ribosomal protein encoding genes from [Fusgr1](https://mycocosm.jgi.doe.gov/Fusgr1/Fusgr1.home.html) and submits BLAST search against all strains.
2. `./submit_codon_optimisation.sh` - submits scripts to estimate various codon usage bias statistics and codon optimisation values.


