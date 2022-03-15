# *Fusarium* Lifestyles

## 8 Selection

1. `qsub gbff_files/ncbi_gbff_download.sh` - downloads GBFF files for the strains used in this study from NCBI; also need [Ilysp1 transcripts downloaded from Mycocosm](https://mycocosm.jgi.doe.gov/Ilysp1/Ilysp1.home.html) in `gbff_files` directory.
2. `./submit_pal2nal.sh` - submits script to pull corresponding nucleotides for all proteins and prepares codon alignments using [PAL2NAL](http://www.bork.embl.de/pal2nal/).
3. `./submit_hyphy.sh` - submits [HyPhy](https://github.com/veg/hyphy) dN/dS methods.
4. :file_folder: `codon_optimisation`
