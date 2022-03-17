# *Fusarium* Lifestyles

## 2 Annotation

### 2.1 Repeatmasking
 
1. `qsub repeatmodeler.sh` - makes custom repeat library for each strain using [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/).
2. `qsub repeatmasker.sh` - uses the custom repeat libraries to softmask assemblies using [RepeatMasker](https://www.repeatmasker.org/RepeatMasker/).
