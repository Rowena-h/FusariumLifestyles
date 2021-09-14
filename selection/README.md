# Fusarium Lifestyles

## 8 Selection

1. `qsub gbff_files/ncbi_gbff_download.sh` - downloads GBFF files from NCBI; also need Ilysp one in `gbff_files` directory.
2. `./submit_pal2nal.sh` - submits script to pull corresponding nucleotides for all proteins and prepares codon alignments.
3. `./submit_hyphy.sh` - submits HyPhy dN/dS selection methods.
4. `codon_optimisation`
