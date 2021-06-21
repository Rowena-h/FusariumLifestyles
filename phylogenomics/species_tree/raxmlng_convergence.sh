module load anaconda3
source activate raxml-ng

raxml-ng --bsconverge --bs-trees fus_proteins_62T.raxml.bootstraps --prefix fus_proteins_62T_convergence_test --seed 2 --threads 1
