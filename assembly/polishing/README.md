# *Fusarium* Lifestyles

## 1 Assembly
### 3 Polishing
 
`qsub polish.sh` - for each strain and assembly tool, maps raw reads to assembly and calculates mapping statistics with BWA-MEM and SAMtools and polishes the assembly with Pilon.

After completing [1 Assembly 4 Assessment](https://github.com/Rowena-h/FusariumLifestyles/tree/main/assembly/assessment) and uploading to NCBI:

`./ncbi_filter.sh` - removes sequences identified as mitochondrial or duplicates by NCBI.
