# *Fusarium* Lifestyles

## 4 Phylogenomics

1. `./submit_alignment.sh` - submits `alignment.sh` for alignment of single copy orthogroups from OrthoFinder with [MAFFT](https://mafft.cbrc.jp/alignment/software/) followed by trimming with [BMGE](https://bmcecolevol.biomedcentral.com/articles/10.1186/1471-2148-10-210) and [trimAl](http://trimal.cgenomics.org/).
2. `./concat.sh` - concatenate single copy orthogroup alignments and prepare partition files using [AMAS](https://github.com/marekborowiec/AMAS).
3. `./submit_modeltestng.sh` - submits `modeltest-ng/modeltestng.sh` to run [ModelTest-NG](https://github.com/ddarriba/modeltest) on all single copy orthogroups (in computationally tractable chunks).
4. `./submit_speciestrees_concatenation.sh` - submits concatenation-based species tree methods - `species_tree/raxmlng.sh` ([RAxML-NG](https://github.com/amkozlov/raxml-ng)) and `species_tree/iqtree.sh` ([IQ-TREE](https://github.com/iqtree/iqtree2)).
5. `./submit_RAxML-NG_genetrees.sh` - submits `RAxMLNG_genetrees.sh` to run RAxML-NG for individual gene trees.
6. `./submit_speciestrees_coalescent.sh` - submits coalescent-based species tree methods - `species_tree/astral.sh` ([ASTRAL-III](https://github.com/smirarab/ASTRAL)) and `species_tree/astral-pro.sh` ([ASTRAL-Pro](https://github.com/chaoszhang/A-pro)) using genes trees.
