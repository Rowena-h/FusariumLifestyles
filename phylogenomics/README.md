# *Fusarium* Lifestyles

## 4 Phylogenomics

1. `./submit_alignment.sh` - submits alignment of single copy orthogroups from OrthoFinder with [MAFFT](https://mafft.cbrc.jp/alignment/software/) followed by trimming with [BMGE](https://bmcecolevol.biomedcentral.com/articles/10.1186/1471-2148-10-210) and [trimAl](http://trimal.cgenomics.org/).
2. `./concat.sh` - concatenates single copy orthogroup alignments and prepares partition files using [AMAS](https://github.com/marekborowiec/AMAS).
3. `./submit_modeltestng.sh` - runs [ModelTest-NG](https://github.com/ddarriba/modeltest) on all single copy orthogroups (in computationally tractable chunks).
4. `./submit_speciestrees_concatenation` - submits concatenation-based species tree methods ([RAxML-NG](https://github.com/amkozlov/raxml-ng) and [IQ-TREE](https://github.com/iqtree/iqtree2)).
5. `./submit_RAxML-NG_genetrees.sh` - submits RAxML-NG for individual gene trees.
6. `./submit_speciestrees_coalescent` - submits coalescent-based species tree methods ([ASTRAL-III](https://github.com/smirarab/ASTRAL) and [ASTRAL-Pro](https://github.com/chaoszhang/A-pro)) using genes trees.