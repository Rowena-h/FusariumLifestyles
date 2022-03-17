# *Fusarium* Lifestyles

## 8 Selection

### 8.2 Codon optimisation

1. `./pull_ribosomes.sh` - extracts ribosomal protein encoding genes from [Fusgr1](https://mycocosm.jgi.doe.gov/Fusgr1/Fusgr1.home.html) and submits `blast.sh` to run BLAST search against all strains in this study.
2. `./submit_codon_optimisation.sh` - submits `codon_optimisation.r` script to estimate various codon usage bias statistics and codon optimisation values.
