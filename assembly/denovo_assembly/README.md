# *Fusarium* Lifestyles

## 1 Assembly

### 1.2 *De novo* genome assembly
 
1. `./submit_assembly.sh` - makes new directory and submits job scripts for each assembly tool - `abyss.sh` ([ABySS](https://github.com/bcgsc/abyss)), `megahit.sh` ([MEGAHIT](https://github.com/voutcn/megahit)) and `spades.sh` ([SPAdes](https://github.com/ablab/spades)).
2. `./abyss_comp.sh` - compares the assembly stats to choose 'best' kmer size for ABySS (must be done after `abyss.sh` has finished for all kmer sizes and strains).
