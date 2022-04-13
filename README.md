# *Fusarium* Lifestyles
 
![Pipeline workflow](pipeline.png)

Bioinformatics analysis pipeline for Hill et al. (2022) Lifestyle transitions in fusarioid fungi are frequent and lack clear genomic signatures. Molecular Biology and Evolution (in press).

The pipeline was written for and run on Queen Mary University of London's [Apocrita HPC facility](http://doi.org/10.5281/zenodo.438045) which uses the Univa Grid Engine batch-queue system. This means that many of the bash scripts (`.sh` file endings) specify core allocation, run times and memory usage allocation that may need to be adapted for different platforms.

---

## 1 Assembly

`cd assembly`

### 1.1 Read quality control

`cd assembly/reads`

Requires raw `fastq.gz` paired-end reads in this directory as well as `TruSeq3-PE.fa` file with adapter sequences downloaded from [here](https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE.fa) (for Illumina NovaSeq 6000 151bp paired-end reads).

1. `qsub trimmomatic.sh` - trims raw reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).
2. `qsub fastqc.sh` - after trimming, checks read quality with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

### 1.2 *De novo* genome assembly

`cd assembly/denovo_assembly`

1. `./submit_assembly.sh` - makes new directory and submits job scripts for each assembly tool - `abyss.sh` ([ABySS](https://github.com/bcgsc/abyss)), `megahit.sh` ([MEGAHIT](https://github.com/voutcn/megahit)) and `spades.sh` ([SPAdes](https://github.com/ablab/spades)).
2. `./abyss_comp.sh` - compares the assembly stats to choose 'best' kmer size for ABySS (must be done after `abyss.sh` has finished for all kmer sizes and strains).

### 1.3 Polishing

`cd assembly/polishing`

1.`qsub polish.sh` - for each strain and assembly tool, maps raw reads to assembly and calculates mapping statistics with [BWA-MEM](https://github.com/lh3/bwa) and [SAMtools](http://www.htslib.org/) and polishes the assembly with [Pilon](https://github.com/broadinstitute/pilon). Also removes sequences <200bp using [Seqtk](https://github.com/lh3/seqtk) for NCBI compliance.

After completing [4 Assessment](https://github.com/Rowena-h/FusariumLifestyles/tree/main/assembly/assessment) and uploading to NCBI:

2.`./ncbi_filter.sh` - removes sequences identified as mitochondrial or duplicates by NCBI (listed in files saved as `duplicates_*` and `mito_remove_*` for each strain) and trims sequences identified as having mitochondrial contaminants (listed in files saved as `mito_trim_*.bed`)  using [bedtools](https://bedtools.readthedocs.io/en/latest/).

### 1.4 Assessment

`cd assembly/assessment`

1. `./submit_assessment.sh` - submits `quast.sh` ([QUAST](https://github.com/ablab/quast)), `busco.sh` ([BUSCO](https://busco.ezlab.org/)) and `blast.sh` ([BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)) scripts for assembly quality statistics. `busco.sh` requires the Hypocreales BUSCO dataset downloaded from [here](https://busco-data.ezlab.org/v4/data/lineages/)).
2. `qsub blobtools.sh` - submits `blobtools.sh` to run [BlobTools](https://github.com/DRL/blobtools) (must be done after `blast.sh` has finished for all strains).

---

## 2 Annotation

`cd annotation`

### 2.1 Repeatmasking

`cd annotation/repeat_masking`

1. `qsub repeatmodeler.sh` - makes custom repeat library for each strain using [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/).
2. `qsub repeatmasker.sh` - uses the custom repeat libraries to softmask assemblies using [RepeatMasker](https://www.repeatmasker.org/RepeatMasker/).

### 2.2 MAKER pipeline

`cd annotation/maker`

Informed by [this](https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2) tutorial. Requires ESTs and proteins from [Fusoxy1](https://mycocosm.jgi.doe.gov/Fusoxy1/Fusoxy1.home.html) and [Fuseq1](https://mycocosm.jgi.doe.gov/Fuseq1/Fuseq1.home.html) ([Mesny et al. 2021](https://doi.org/10.1038/s41467-021-27479-y)) downloaded from Mycocosm in this directory.

#### Round 1

`cd annotation/maker/round1`

`qsub maker.sh` - submits first run of [MAKER](http://www.yandell-lab.org/software/maker.html) using ESTs and proteins, as indicated in `.ctl` files.

#### Round 2

`cd annotation/maker/round2`

1. `qsub training_snap/snap.sh` - trains [SNAP](https://github.com/KorfLab/SNAP) using gene models from the first MAKER round.
2. `qsub maker2.sh` - submits second run of MAKER using trained SNAP (as indicated in `.ctl` files).

#### Round 3

`cd annotation/maker/round3`

1. `qsub training_snap2/snap2.sh` - trains SNAP again using gene models from the second MAKER round.
2. `qsub maker3.sh` - submits third run of MAKER using second trained SNAP (as indicated in `.ctl` files).
3. `qsub rename.sh` - after obtaining unique locus tags from e.g. NCBI (see `locus_tags.txt`), renames IDs in gff and fasta files.
4. `qsub gag.sh` - runs [GAG](https://github.com/genomeannotation/GAG/) to remove introns <10bp, remove terminal Ns and correct start and stop codons in gff file for NCBI compliance.

---

## 3 Orthology inference

`cd orthology_inference`

1. `./submit_protein_download.sh` - submits `ncbi_ftp_links.r` and `protein_download.sh` to download of predicted protein sets of *Fusarium* strains from NCBI.
2. `qsub orthofinder.sh` - submits orthology inference using [OrthoFinder](https://github.com/davidemms/OrthoFinder).

---	

## 4 Phylogenomics

`cd phylogenomics`

1. `./submit_alignment.sh` - submits `alignment.sh` for alignment of single copy orthogroups from OrthoFinder with [MAFFT](https://mafft.cbrc.jp/alignment/software/) followed by trimming with [BMGE](https://bmcecolevol.biomedcentral.com/articles/10.1186/1471-2148-10-210) and [trimAl](http://trimal.cgenomics.org/).
2. `./concat.sh` - concatenate single copy orthogroup alignments and prepare partition files using [AMAS](https://github.com/marekborowiec/AMAS).
3. `./submit_modeltestng.sh` - submits `modeltest-ng/modeltestng.sh` to run [ModelTest-NG](https://github.com/ddarriba/modeltest) on all single copy orthogroups (in computationally tractable chunks).
4. `./submit_speciestrees_concatenation.sh` - submits concatenation-based species tree methods - `species_tree/raxmlng.sh` ([RAxML-NG](https://github.com/amkozlov/raxml-ng)) and `species_tree/iqtree.sh` ([IQ-TREE](https://github.com/iqtree/iqtree2)).
5. `./submit_RAxML-NG_genetrees.sh` - submits `RAxMLNG_genetrees.sh` to run RAxML-NG for individual gene trees.
6. `./submit_speciestrees_coalescent.sh` - submits coalescent-based species tree methods - `species_tree/astral.sh` ([ASTRAL-III](https://github.com/smirarab/ASTRAL)) and `species_tree/astral-pro.sh` ([ASTRAL-Pro](https://github.com/chaoszhang/A-pro)) using genes trees.

---

## 5 Divergence time estimation

`cd divergence_time_estimation`

### 5.1 SortaDate

`cd divergence_time_estimation/sortadata`

`qsub sortadate` - reroots gene and RAxML-NG species tree and runs with [SortaDate](https://github.com/FePhyFoFum/SortaDate) to filter for top ten 'clock-like' genes.

### 5.2 MCMCTree

`cd divergence_time_estimation/mcmctree`

1. `qsub mcmctree_dating_step1.sh` - adds secondary time calibrations to species tree nodes and submits first step of approximate likelihood divergence time estimation with protein data using MCMCTree from [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) (see [tutorial](http://abacus.gene.ucl.ac.uk/software/MCMCtree.Tutorials.pdf)).
2. `Rscript estimate_rate.r` - estimates the scaling parameter for the substitution rate prior to be added to `mcmctree_step2_independent.ctl` and `mcmctree_step2_correlated.ctl`.
3. `./submit_mcmctree_dating_step2.sh` - submits `mcmctree_independent.sh` and `mcmctree_correlated.sh` for second step of approximate likelihood estimation for both independent and correlated rates relaxed clock models.

---

## 6 CSEP & CAZyme prediction

`cd CSEP_CAZyme_prediction`

1. `./submit_CSEPprediction.sh` - submits all programmes in the CSEP prediction pipeline - `signalp/signalp.sh` ([SignalP](https://services.healthtech.dtu.dk/service.php?SignalP-5.0)), `targetp/targetp.sh` ([TargetP](https://services.healthtech.dtu.dk/service.php?TargetP-2.0)), `phobius/phobius.sh` ([Phobius](https://phobius.sbc.su.se/instructions.html)), `tmhmm/tmhmm.sh` ([TMHMM](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0)), `prosite/ps_scan.sh` ([ps_scan](https://prosite.expasy.org/scanprosite/)), `nucpred/nucpred.sh` ([NucPred](https://nucpred.bioinfo.se/nucpred/)), `predgpi/predgpi.sh` which in turn submits `predgpi/PredGPI.r` to use the R package [ragp](https://rdrr.io/github/missuse/ragp/man/get_pred_gpi.html) ([PredGPI](http://gpcr.biocomp.unibo.it/predgpi/)) and  `effectorp/effectorp.sh` ([EffectorP](https://github.com/JanaSperschneider/EffectorP-3.0)).
2. `./submit_CSEPfilter.sh` - submits `CSEPfilter` to produce lists of CSEPs from all programme results.
3. `./submit_CSEPblast.sh` -  submits `blastp/blastp.sh` to BLAST of CSEPs against the [PHI-base database](http://www.phi-base.org/) (requires `phi-base_current.csv` and `phi-base_current.fas` to be downloaded from [here](http://www.phi-base.org/downloadLink.htm) into the `blastp` directory).
4. `./submit_CAZymeprediction.sh` - submits `run_dbcan.sh/run_dbcan.sh` to run run_dbcan. 
5. `qsub submit_orthogroupparsing.sh` - submits `orthogroup_parser.r` to make abundance matrices of orthogroups for all strains and categorises whether they are CSEPs/CAZymes and core/accessory/specific.

---

## 7 Lifestyle comparison

`cd lifestyle_comparison`

`./submit_lifestyletest.sh` - submits `lifestyle_v_phylogeny.r` to prepare input file for `lifestyle-test.sh` which runs PERMANOVA-based lifestyle test on orthogroup and CSEP presence absence matrices; `run_edited.py` is modified from the original script `run.py` by [Mesny & Vannier](https://github.com/fantin-mesny/Effect-Of-Biological-Categories-On-Genomes-Composition).

---

## 8 Selection

`cd selection`

### 8.1 dN/dS analysis

1. `qsub gbff_files/ncbi_gbff_download.sh` - downloads GBFF files for the strains used in this study from NCBI; also need [Ilysp1 transcripts downloaded from Mycocosm](https://mycocosm.jgi.doe.gov/Ilysp1/Ilysp1.home.html) in `gbff_files` directory.
2. `./submit_pal2nal.sh` - submits `pal2nal.sh` script to pull corresponding nucleotides for all proteins using `pull_nucleotides.py` and prepares codon alignments using [PAL2NAL](http://www.bork.embl.de/pal2nal/).
3. `./submit_hyphy.sh` - prepares file inputs and submits scripts for [HyPhy](https://github.com/veg/hyphy) dN/dS methods - `hyphy/BUSTED.sh`, `hyphy/aBSREL.sh` and `hyphy/Contrast-FEL.sh`.

### 8.2 Codon optimisation

`cd selection/codon_optimisation`

1. `./pull_ribosomes.sh` - extracts ribosomal protein encoding genes from [Fusgr1](https://mycocosm.jgi.doe.gov/Fusgr1/Fusgr1.home.html) and submits `blast.sh` to run BLAST search against all strains in this study.
2. `./submit_codon_optimisation.sh` - submits `codon_optimisation.r` script to estimate various codon usage bias statistics and codon optimisation values.

---

## 9 Statistics and data visualisation

`Rscript stats_and_plots.r`
