# *Fusarium* Lifestyles

## 5 Divergence time estimation

### 5.2 MCMCTree

1. `qsub mcmctree_dating_step1.sh` - adds secondary time calibrations to species tree nodes and submits first step of approximate likelihood divergence time estimation with protein data using MCMCTree from [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) (see [tutorial](http://abacus.gene.ucl.ac.uk/software/MCMCtree.Tutorials.pdf)).
2. `Rscript estimate_rate.r` - estimates the scaling parameter for the substitution rate prior to be added to `mcmctree_step2_independent.ctl` and `mcmctree_step2_correlated.ctl`.
3. `./submit_mcmctree_dating_step2.sh` - submits `mcmctree_independent.sh` and `mcmctree_correlated.sh` for second step of approximate likelihood estimation for both independent and correlated rates relaxed clock models.
