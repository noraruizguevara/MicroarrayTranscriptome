# MicroarrayTranscriptome
Codes for transcriptome-microarray analysis

# PASO 0: ¿De qué trata este repositorio?
```r
En este repositorio se presentan comandos para la identificacion de los genes y proteinas que forman parte de la libreria
de microarreglos empleada en el "analisis previo".

```

# PASO 1: BLAST vs NR-NCBI, obtencion de archivos (csv y tsv)
```r
## Con el objetivo de identificar los "contigs" obtenidos con el ensayo microarray
## se procede a realizar un analisis BLASTx con tres formatos distintos

blastx -query CONTIG.fasta -db NR -out CONTIG_BLASTX2nr -outfmt 5 -evalue 0.0001 -gapopen 11 -gapextend 1 -word_size 3 -matrix BLOSUM62 -max_target_seqs 20 -num_threads 30

blastx -query CONTIG.fasta -db NR -out CONTIG_csvBLASTX2nr -outfmt 10 -evalue 0.0001 -gapopen 11 -gapextend 1 -word_size 3 -matrix BLOSUM62 -max_target_seqs 20 -num_threads 30

blastx -query CONTIG.fasta -db NR -out CONTIG_txtBLASTX2nr -outfmt 0 -evalue 0.0001 -gapopen 11 -gapextend 1 -word_size 3 -matrix BLOSUM62 -num_descriptions 20 -num_alignments 20 -num_threads 30
```

# PASO 2: Identificar los mejores HITS (proteinas) para cada uno de los transcritos
```r
## El archivo "CONTIG_csvBLASTX2nr" contiene informacion para identificar
## los mejores HITS (aquellos con mayor %identidad y %cobertura)
## las accesiones correspondientes son utiles para extraer sus respectivos
## peptidos a traves de "batch entrez" (https://www.ncbi.nlm.nih.gov/sites/batchentrez)

# 1. setear el directorio de trabajo
setwd("C:/Users/HP/Documents/lepidium/BLAST/contigs_vs_NR/uncompress")
dir()

# 2. leer "CONTIG_csvBLASTX2nr"
data <- read.csv("CONTIG_csvBLASTX2nr", header=F)
head(data)

# 3. añadir headers
names(data) <- c("qseqid","sseqid","pident","length",
"mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
head(data)
dim(data)

# 4. leer el archivo de longitud (nt) de cada secuencia de "CONTIG.fasta"
#    transformarlo en un data.frame de dos columnas con nombres
length <- read.csv("lengths.txt", header=F)
dim(length)
head(length)
names(length) <- c("qseqid","qbases")

# 5. identificar todos los hits unicos y almacenarlos en un
#    archivo de nombre "hits.txt"

hist(data$pident)
hist(data$length)
View(data)

length(sort(unique(data$sseqid)))
hits <- sort(unique(data$sseqid))
write.table(hits,"hits.txt",col.names=F,row.names=F,quote=F)

# como hay 20 hits por secuencia, entonces: 
7158*20

# 6. Pero nosotros solo necesitamos el mejor hit (1)
#    para ello nos basaremos en el parametro "bit score"
#    ordenamos los resultados de cada transcrito "subject"
#    por la columna "bit score" y solo retemos el que
#    corresponde al mayor valor.

length(sort(unique(data$qseqid)))

a1 <- c() ; 
a2 <- c() ; 
a3 <- c() ; 
a4 <- c() ; 

contig <- unique(data$qseqid)
for (i in contig){
a1 <- data[data$qseqid %in% i, ]
a2 <- a1[order(-a1$bitscore),]
a3 <- a2[1,]
a4 <- rbind(a4,a3)
}

dim(a4)
a4[1:20,]

# 7. cambiamos los nombres de las filas y graficamos
#    histogramas exploratorios, luego generamos un nuevo
#    output de hits "hits_2.txt"

row_selected <- row.names(a4)

hist(a4$pident)
hist(a4$length)
hist(a4$mismatch)
hist(a4$evalue)

hits <- sort(unique(a4$sseqid))
length(sort(unique(a4$sseqid)))

write.table(hits,"hits_2.txt",col.names=F,row.names=F,quote=F)

# 8. Una verificacion visual del resultado identifica accesiones
#    con nombres no identificables, seran almacenadas en el
#    objeto "exc".  

exc <- c("pir|B84683|","pir|C84609|","pir|T45670|","prf||1804333D")

head(data)

d1 <- unique(data[data$sseqid %in% exc, "qseqid"])

a1 <- c() ; 
a2 <- c() ; 
a3 <- c() ; 
a7 <- c() ; 

contig <- d1
for (i in contig){
a1 <- data[data$qseqid %in% i, ]
a2 <- a1[order(-a1$bitscore),]
a3 <- a2[1:3,]
a7 <- rbind(a7,a3)
}

dim(a7)
a7

# 9. Las siguientes lineas tienen el objetivo de identificar y corregir manualmente
#    los hits con codigos de accesion ilegibles, el resultado se almacenará
#    en el objeto "best.hit.per.contig.csv"

todas.filas <- row_selected
excluir.filas <- c("14690","17919","117572","120611")
incluir.filas <- c("14691","17915","117573","120612")
intersect(todas.filas,excluir.filas)

j1 <- setdiff(todas.filas,excluir.filas)
length(j1)
j2 <- sort(c(j1,incluir.filas))
length(j2)

j3 <- as.character(sort(as.numeric(j2)))

data2 <- data[j3,]
dim(data2)
head(data2)

write.table(data2,"best.hit.per.contig.csv", row.names=F,quote=F,sep=",")

# 10. generar un output con los nuevos hits recuperados

hits.new <- sort(unique(data2$sseqid))
nuevos <- data[incluir.filas,"sseqid"]
write.table(nuevos,"hits_4.txt",col.names=F,row.names=F,quote=F)

## Los BEST HITS se pueden descargar y anotar, ello nos permitirá
## identificar las caracteristicas de los transcritos
```

# PASO 2: Anotación de los contigs con TRANSDECODER + EGGNOGMAPPER 
```r
## Otra forma de obtener informacion de los transcritos
## es inferir los ORFs, péptidos y anotarlos con EGGNOG MAPPER 

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
