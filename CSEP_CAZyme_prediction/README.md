# *Fusarium* Lifestyles

## 6 CSEP prediction

1. `./submit_CSEPprediction.sh` - submits all programmes in the CSEP prediction pipeline - [SignalP](https://services.healthtech.dtu.dk/service.php?SignalP-5.0), [TargetP](https://services.healthtech.dtu.dk/service.php?TargetP-2.0), [Phobius](https://phobius.sbc.su.se/instructions.html), [TMHMM](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0), [ps_scan](https://prosite.expasy.org/scanprosite/), [NucPred](https://nucpred.bioinfo.se/nucpred/), [PredGPI](http://gpcr.biocomp.unibo.it/predgpi/) (accessed via the R package [ragp](https://rdrr.io/github/missuse/ragp/man/get_pred_gpi.html)) and [EffectorP](https://github.com/JanaSperschneider/EffectorP-3.0).
2. `./submit_CSEPfilter.sh` - runs all programme results through `CSEPfilter` to produce lists of CSEPs.
3. `./submit_CSEPblast.sh` -  submits BLAST of CSEPs against the [PHI-base database](http://www.phi-base.org/) (requires `phi-base_current.csv` and `phi-base_current.fas` to be downloaded from [here](http://www.phi-base.org/downloadLink.htm) into the `blastp` directory).
4. `./submit_CAZymeprediction.sh` - submits run_dbcan for all strains in this study. 
5. `qsub submit_orthogroupparsing.sh` - makes abundance matrices of orthogroups for all strains and categorises whether they are CSEPs/CAZymes and core/accessory/specific.