# Fusarium Lifestyles

## 2 Annotation
### 2 MAKER pipeline
#### Round 3

1. `qsub training_snap2/snap2.sh` - trains SNAP again using gene models from the second MAKER round.
2. `qsub maker3.sh` - third run of MAKER using second trained SNAP (as indicated in `.ctl` files).
3. `qsub rename.sh` - after obtaining unique locus tags from e.g. NCBI, renames IDs in gff and fasta files.
4. `qsub gag.sh` - removes introns <10bp, removes terminal Ns and corrects start and stop codons in gff file for NCBI compliance.
