# *Fusarium* Lifestyles

## 1 Assembly

### 1.1 Read quality control
 
Requires raw `fastq.gz` paired-end reads in this directory as well as `TruSeq3-PE.fa` file with adapter sequences downloaded from [here](https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE.fa) (for Illumina NovaSeq 6000 151bp paired-end reads).

1. `qsub trimmomatic.sh` - trims raw reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).
2. `qsub fastqc.sh` - after trimming, checks read quality with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).