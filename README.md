# MicroarrayTranscriptome
Codes for transcriptome-microarray analysis

# PASO 0: ¿De qué trata este repositorio?
```r
```

# PASO 1: BLAST vs NR-NCBI 
```r
Con el objetivo de identificar los contigs obtenidos con el ensayo microarray
se procede a realizar un analisis BLASTx con tres formatos distintos

blastx -query CONTIG.fasta -db NR -out CONTIG_BLASTX2nr -outfmt 5 -evalue 0.0001 -gapopen 11 -gapextend 1 -word_size 3 -matrix BLOSUM62 -max_target_seqs 20 -num_threads 30

blastx -query CONTIG.fasta -db NR -out CONTIG_csvBLASTX2nr -outfmt 10 -evalue 0.0001 -gapopen 11 -gapextend 1 -word_size 3 -matrix BLOSUM62 -max_target_seqs 20 -num_threads 30

blastx -query CONTIG.fasta -db NR -out CONTIG_txtBLASTX2nr -outfmt 0 -evalue 0.0001 -gapopen 11 -gapextend 1 -word_size 3 -matrix BLOSUM62 -num_descriptions 20 -num_alignments 20 -num_threads 30

```
