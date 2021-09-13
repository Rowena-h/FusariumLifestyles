# Fusarium Lifestyles

## 2 Annotation
### 2 MAKER pipeline
#### Round 3

1. `qsub training_snap2/snap2.sh` - train SNAP again using gene models from the second MAKER round.
2. `qsub maker3.sh` - third run of MAKER using second trained SNAP (as indicated in `.ctl` files).

