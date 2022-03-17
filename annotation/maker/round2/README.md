# *Fusarium* Lifestyles

## 2 Annotation

### 2.2 MAKER pipeline

#### Round 2

1. `qsub training_snap/snap.sh` - trains [SNAP](https://github.com/KorfLab/SNAP) using gene models from the first MAKER round.
2. `qsub maker2.sh` - submits second run of MAKER using trained SNAP (as indicated in `.ctl` files).
