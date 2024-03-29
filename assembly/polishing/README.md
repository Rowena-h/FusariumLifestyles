# *Fusarium* Lifestyles

## 1 Assembly

### 1.3 Polishing
 
1.`qsub polish.sh` - for each strain and assembly tool, maps raw reads to assembly and calculates mapping statistics with [BWA-MEM](https://github.com/lh3/bwa) and [SAMtools](http://www.htslib.org/) and polishes the assembly with [Pilon](https://github.com/broadinstitute/pilon). Also removes sequences <200bp using [Seqtk](https://github.com/lh3/seqtk) for NCBI compliance.

After completing [4 Assessment](https://github.com/Rowena-h/FusariumLifestyles/tree/main/assembly/assessment) and uploading to NCBI:

2.`./ncbi_filter.sh` - removes sequences identified as mitochondrial or duplicates by NCBI (listed in files saved as `duplicates_*` and `mito_remove_*` for each strain) and trims sequences identified as having mitochondrial contaminants (listed in files saved as `mito_trim_*.bed`)  using [bedtools](https://bedtools.readthedocs.io/en/latest/).
