# *Fusarium* Lifestyles

## 1 Assembly

### 1.4 Assessment
 
1. `./submit_assessment.sh` - submits `quast.sh` ([QUAST](https://github.com/ablab/quast)), `busco.sh` ([BUSCO](https://busco.ezlab.org/)) and `blast.sh` ([BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)) scripts for assembly quality statistics. `busco.sh` requires the Hypocreales BUSCO dataset downloaded from [here](https://busco-data.ezlab.org/v4/data/lineages/)).
2. `qsub blobtools.sh` - submits `blobtools.sh` to run [BlobTools](https://github.com/DRL/blobtools) (must be done after `blast.sh` has finished for all strains).
