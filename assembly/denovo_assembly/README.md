# Fusarium Lifestyles

## 1 Assembly
### 2 *De novo* genome assembly
 
1. `./submit_assembly` - makes new directory and submits job scripts for each assembly tool (ABySS, MEGAHIT, SPAdes).
2. `./abyss_comp.sh` - after ABySS has finished running for all kmer sizes, compare the assembly stats to choose 'best' kmer size.