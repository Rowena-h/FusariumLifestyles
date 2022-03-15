# *Fusarium* Lifestyles

## 2 Annotation
### 1 Repeatmasking
 
1. `qsub repeatmodeler.sh` - makes custom repeat library for each strain using [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/).
2. `qsub repeatmasker.sh` - after repeat modelling, uses custom repeat library to softmask the assembly using [RepeatMasker](https://www.repeatmasker.org/RepeatMasker/).
