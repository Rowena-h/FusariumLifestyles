#!/bin/sh

#Search annotation for ribosomal proteins
grep "ribosomal protein" Fusgr1_GeneCatalog_proteins_20110524_IPR.tab > ribosomal_proteins.tsv
#Remove mitochondrial ones
sed -i '/Mitochondrial/d' ribosomal_proteins.tsv	

#Format ribosomal protein file for bedtools
cut -f1,8,9 ribosomal_proteins.tsv > ribosomal_proteins.bed
sed -i 's/,//g' ribosomal_proteins.bed
sed -i -e 's/jgi|Fusgr1|//' Fusgr1_GeneCatalog_proteins_20110524.aa.fasta
sed -i -e 's/|.*//' Fusgr1_GeneCatalog_proteins_20110524.aa.fasta

module load bedtools

#Pull sequences for ribosomal proteins
bedtools getfasta -fi Fusgr1_GeneCatalog_proteins_20110524.aa.fasta -bed ribosomal_proteins.bed -fo ribosomal_proteins.fasta

#Format fasta for blasting
sed -i 's/:.*//' ribosomal_proteins.fasta
awk -F '\t' '{print $1 "\t" $3}' ribosomal_proteins.tsv | sed "s/[^a-zA-Z0-9]/_/g" | sed 's/_/\t/' > headers
awk -F '\t' 'FNR==NR {f2[$1]=$2;next} $2 in f2 {$2=f2[$2]}1' headers FS='>' OFS='>' ribosomal_proteins.fasta > tmp && mv tmp ribosomal_proteins.fasta

rm headers

NUM_TAXA=$(ls ../../proteins/*.faa | sed 's#\.\./\.\./proteins/##' | wc -l)

qsub -t 1-${NUM_TAXA} blast.sh
