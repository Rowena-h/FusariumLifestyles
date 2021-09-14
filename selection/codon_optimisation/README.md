# Fusarium Lifestyles

## 8 Selection
### Codon optimisation

1. `./pull_ribosomes.sh` - extracts ribosomal protein encoding genes from [Fusgr1](https://mycocosm.jgi.doe.gov/Fusgr1/Fusgr1.home.html) and submits BLAST search against all strains.
2. `./submit_codon_optimisation.sh` - submits scripts to estimate various codon usage bias statistics and codon optimisation values.
