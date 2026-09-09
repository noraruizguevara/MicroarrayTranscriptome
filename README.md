# MicroarrayTranscriptome
Codes for transcriptome-microarray analysis

# PASO 0: ¿De qué trata este repositorio?
```r
En este repositorio se presentan comandos para la identificacion de los genes y proteinas que forman parte de la libreria
de microarreglos empleada en el "analisis previo".

```

# PASO 1: BLAST vs NR-NCBI 
```r
## Con el objetivo de identificar los "contigs" obtenidos con el ensayo microarray
## se procede a realizar un analisis BLASTx con tres formatos distintos

blastx -query CONTIG.fasta -db NR -out CONTIG_BLASTX2nr -outfmt 5 -evalue 0.0001 -gapopen 11 -gapextend 1 -word_size 3 -matrix BLOSUM62 -max_target_seqs 20 -num_threads 30

blastx -query CONTIG.fasta -db NR -out CONTIG_csvBLASTX2nr -outfmt 10 -evalue 0.0001 -gapopen 11 -gapextend 1 -word_size 3 -matrix BLOSUM62 -max_target_seqs 20 -num_threads 30

blastx -query CONTIG.fasta -db NR -out CONTIG_txtBLASTX2nr -outfmt 0 -evalue 0.0001 -gapopen 11 -gapextend 1 -word_size 3 -matrix BLOSUM62 -num_descriptions 20 -num_alignments 20 -num_threads 30
```

# PASO 2: 
```r
```

# PASO 2: Anotación de los contigs con TRANSDECODER + EGGNOGMAPPER 
```r

#!/usr/bin/bash

## 1. instalar TRANSDECODER ##
conda create -n transdecoder ;
conda activate transdecoder ;
conda install bioconda::transdecoder ;

## 2. identificar peptidos mayores o iguales a 30aa ##
TransDecoder.LongOrfs -t contigs.array.fasta -m 30 ;

## 3. filtrar los datos anteriores y retener los peptidos de mejor calidad ##
TransDecoder.Predict -t contigs.array.fasta -T 1000 ;

## 4. chequear resultados ## 
grep "^>" contigs.array.fasta.transdecoder.pep | wc -l ;
 # 9066
grep "^>" contigs.array.fasta | wc -l ;
# 7158

## 5. anotar con EGGNOG-MAPPER ##
emapper.py -i contigs.array.fasta.transdecoder.clean.pep \
           -o contigs_eggnog_annotation \
           --cpu 28 \
           -m diamond \
           --dmnd_db /mnt/d/TESIS_MACA_2026/BIOINFORMATICS/4.EGGNOG.ANNOTATION/eggnog_proteins.dmnd \
           --data_dir /mnt/d/TESIS_MACA_2026/BIOINFORMATICS/4.EGGNOG.ANNOTATION/ \
           --sensmode fast ;
ls ;
```
