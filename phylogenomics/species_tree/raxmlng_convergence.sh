module load anaconda3
source activate raxml-ng

raxml-ng --bsconverge --bs-trees raxml-ng/fus_proteins_62T.raxml.bootstraps --prefix raxml-ng/fus_proteins_62T_convergence_test --seed 2 --threads 1
