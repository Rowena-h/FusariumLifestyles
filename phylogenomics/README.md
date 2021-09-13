# Fusarium Lifestyles

## 4 Phylogenomics

1. `./submit_alignment.sh` - submits alignment of single copy orthogroups from OrthoFinder with MAFFT followed by trimming with BMGE.
2. `./concat.sh` - concatenate single copy orthogroup alignments and prepare partition files.
3. `./submit_modeltestng.sh` - runs ModelTest-NG on all single copy orthogroups (in computationally tractable chunks).
4. `./submit_speciestrees_concatenation` - submits concatenation-based species tree methods (RAxML-NG and IQ-TREE).
5. `./submit_RAxML-NG_genetrees.sh` - submits RAxML-NG for individual gene trees.
6. `qsub species_tree/astral.sh` - submits coalescent-based species tree method (ASTRAL-III) using genes trees.