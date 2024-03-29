# *Fusarium* Lifestyles

## 6 CSEP & CAZyme prediction

1. `./submit_CSEPprediction.sh` - submits all programmes in the CSEP prediction pipeline - `signalp/signalp.sh` ([SignalP](https://services.healthtech.dtu.dk/service.php?SignalP-5.0)), `targetp/targetp.sh` ([TargetP](https://services.healthtech.dtu.dk/service.php?TargetP-2.0)), `phobius/phobius.sh` ([Phobius](https://phobius.sbc.su.se/instructions.html)), `tmhmm/tmhmm.sh` ([TMHMM](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0)), `prosite/ps_scan.sh` ([ps_scan](https://prosite.expasy.org/scanprosite/)), `nucpred/nucpred.sh` ([NucPred](https://nucpred.bioinfo.se/nucpred/)), `predgpi/predgpi.sh` which in turn submits `predgpi/PredGPI.r` to use the R package [ragp](https://rdrr.io/github/missuse/ragp/man/get_pred_gpi.html) ([PredGPI](http://gpcr.biocomp.unibo.it/predgpi/)) and  `effectorp/effectorp.sh` ([EffectorP](https://github.com/JanaSperschneider/EffectorP-3.0)).
2. `./submit_CSEPfilter.sh` - submits `CSEPfilter` to produce lists of CSEPs from all programme results.
3. `./submit_CSEPblast.sh` -  submits `blastp/blastp.sh` to BLAST of CSEPs against the [PHI-base database](http://www.phi-base.org/) (requires `phi-base_current.csv` and `phi-base_current.fas` to be downloaded from [here](http://www.phi-base.org/downloadLink.htm) into the `blastp` directory).
4. `./submit_CAZymeprediction.sh` - submits `run_dbcan/run_dbcan.sh` to run [run_dbcan](https://github.com/linnabrown/run_dbcan). 
5. `qsub submit_orthogroupparsing.sh` - submits `orthogroup_parser.r` to make abundance matrices of orthogroups for all strains and categorises whether they are CSEPs/CAZymes and core/accessory/specific.
