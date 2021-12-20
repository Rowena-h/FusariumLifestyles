# *Fusarium* Lifestyles

## 5 Divergence time estimation
### 2 MCMCTree

1. `qsub submit_mcmctree_dating_step1.sh` - submits first step of approximate likelihood divergence time estimation with protein data (see [tutorial](http://abacus.gene.ucl.ac.uk/software/MCMCtree.Tutorials.pdf).
2. `Rscript estimate_rate.r` - estimates the scaling parameter for the substitution rate prior using the species tree.
3. `./submit_mcmctree_dating_step2.sh` - submitssecond step of approximate likelihood estimation for both independent and correlated rates relaxed clock models.
