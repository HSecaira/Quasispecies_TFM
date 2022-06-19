#!/bin/bash


# Mamba es similar a Anaconda pero est ́a implementada en C++
#conda install -y mamba

# Crear un ambiente con cualquier nombre e instalar las herramientas
#mamba create -y -n SeqsAnalyses cutadapt fastqc multiqc \
#bwa samtools igv qualimap flash2 parallel fastx_toolkit

# Activar el ambiente
source activate SeqsAnalyses

# Directorio de trabajo
WD=~/SeqsProcess
mkdir -p $WD

# Y otros directorios
mkdir -p $WD/data/Reg1 $WD/res/Reg1

# Crear directorios
mkdir -p $WD/out/fastqc/Reg1/before_trimming $WD/out/multiqc/Reg1/before_trimming  
mkdir -p $WD/log/fastqc/Reg1/before_trimming $WD/log/multiqc/Reg1/before_trimming

# Generar un reporte individual para cada archivo con FastQC - R1
find $WD/res/Reg1 -name "*R1*" | sed 's/.fastq.gz$//' | parallel -j $(nproc) \
	"fastqc $WD/res/Reg1/{/.}.fastq.gz -o \
	$WD/out/fastqc/Reg1/before_trimming 2>&1 | tee \
	$WD/log/fastqc/Reg1/before_trimming/{/.}.log"

# Generar un reporte individual para cada archivo con FastQC - R2
find $WD/res/Reg1 -name "*R2*" | sed 's/.fastq.gz$//' | parallel -j $(nproc) \
	"fastqc $WD/res/Reg1/{/.}.fastq.gz -o \
	$WD/out/fastqc/Reg1/before_trimming 2>&1 | tee \
	$WD/log/fastqc/Reg1/before_trimming/{/.}.log"

# Generar un solo reporte con MultiQC 
multiqc "$WD/out/fastqc/Reg1/before_trimming/" -f -m fastqc -o \
	"$WD/out/multiqc/Reg1/before_trimming" 2>&1 | tee \
	"$WD/log/multiqc/Reg1/before_trimming/Report.log"

# Crear directorios
mkdir -p $WD/out/cutadapt/Reg1 $WD/out/cutadapt/too_short/Reg1 \
	$WD/out/cutadapt/untrimmed/Reg1 $WD/log/cutadapt/Reg1

# Secuencias de primers
PR1=GAATGTTGGTGACATACTTGCT
PR2=TTGCCATGATCAAATTGACC
# Longitud de barcode
LB=10

find $WD/res/Reg1 -name "*R1_001.fastq.gz" | sed 's/_R1_001.fastq.gz$//' | parallel \
	-j $(nproc) "cutadapt -q 25 -g ^$PR1 -O 13 \
	-G ^N\{$LB\}$PR2 -m 50 -o $WD/out/cutadapt/Reg1/{/.}_R1_001_all_trim.fastq.gz \
	-p $WD/out/cutadapt/Reg1/{/.}_R2_001_all_trim.fastq.gz \
	--pair-filter=any --too-short-output \
	$WD/out/cutadapt/too_short/Reg1/{/.}_R1_001_less_m50.fastq \
	--too-short-paired-output \
	$WD/out/cutadapt/too_short/Reg1/{/.}_R2_001_less_m50.fastq \
	--untrimmed-output \
	$WD/out/cutadapt/untrimmed/Reg1/{/.}_R1_001_untrimmed.fastq \
	--untrimmed-paired-output \
	$WD/out/cutadapt/untrimmed/Reg1/{/.}_R2_001_untrimmed.fastq \
	$WD/res/Reg1/{/.}_R1_001.fastq.gz $WD/res/Reg1/{/.}_R2_001.fastq.gz  \
	2>&1 | tee $WD/log/cutadapt/Reg1/{/.}_001.log"

# Crear directorios
mkdir -p $WD/log/multiqc/cutadapt/Reg1 $WD/out/multiqc/Cutadapt/Reg1 

# Reporte de calidad
multiqc "$WD/log/cutadapt/Reg1/" -f -v -m cutadapt -o \
	"$WD/out/multiqc/cutadapt/Reg1" 2>&1 | tee \
	"$WD/log/multiqc/cutadapt/Reg1/Report.log"

# Crear directorios
mkdir -p $WD/out/fastqc/Reg1/after_trimming $WD/out/multiqc/Reg1/after_trimming \
	$WD/log/multiQC/Reg1/after_trimming $WD/log/fastqc/Reg1/after_trimming

# Reporte de calidad con fastqc
find $WD/out/cutadapt/Reg1 -name "*.fastq.gz" | sed 's/.fastq.gz$//' | parallel \
	-j $(nproc) "fastqc $WD/out/cutadapt/Reg1/{/.}.fastq.gz \
	-o $WD/out/fastqc/Reg1/after_trimming 2>&1 | tee \
	$WD/log/fastqc/Reg1/after_trimming/{/.}.log"

# Reporte de calidad con multiqc
multiqc "$WD/out/fastqc/Reg1/after_trimming/" -f -m fastqc -o \
	"$WD/out/multiqc/Reg1/after_trimming " 2>&1 | tee \
	"$WD/log/multiQC/Reg1/after_trimming/Report.log"

# Crear directorios
mkdir -p $WD/out/flash2/Reg1 $WD/log/flash2/Reg1

# Unir lecturas
find $WD/out/cutadapt/Reg1 -name "*_R1_001_all_trim.fastq.gz" | sed \
	's/_R1_001_all_trim.fastq.gz$//' | parallel -j $(nproc) "flash2 -z -o \
	{/.}_001_all_trim_merged -d $WD/out/flash2/Reg1 -M 180 \
	$WD/out/cutadapt/Reg1/{/.}_R1_001_all_trim.fastq.gz \
	$WD/out/cutadapt/Reg1/{/.}_R2_001_all_trim.fastq.gz 2>&1 | tee \
	$WD/log/flash2/Reg1/{/.}_001_all_trim.log"

# Crear directorios
mkdir -p $WD/log/multiqc/flash2/Reg1 $WD/out/multiqc/flash2/Reg1

# Mover los archivos .hist al directorio log
for f in "$WD/out/flash2/Reg1/*.hist"; do mv $f "$WD/log/flash2/Reg1/"; done;

# Reporte de calidad
multiqc "$WD/log/flash2/Reg1/" -f -m flash -o \
	"$WD/out/multiqc/flash2/Reg1" 2>&1 | tee \
	"$WD/log/multiqc/flash2/Reg1/Report.log"

# Crear directorios
mkdir -p $WD/out/fastqc/Reg1/after_merging $WD/log/fastqc/Reg1/after_merging \
	$WD/out/multiqc/Reg1/after_merging $WD/log/multiqc/Reg1/after_merging

# Reporte de calidad con fastqc
find $WD/out/flash2/Reg1/ -name "*.extendedFrags.fastq.gz" | sed \
	's/.extendedFrags.fastq.gz$//' | parallel -j $(nproc) \
	"fastqc $WD/out/flash2/Reg1/{/.}.extendedFrags.fastq.gz -o \
	$WD/out/fastqc/Reg1/after_merging 2>&1 | tee $WD/log/fastqc/Reg1/after_merging/{/.}.log"

# Reporte de calidad con multiqc
multiqc "$WD/out/fastqc/Reg1/after_merging/" -f -m fastqc -o \
	"$WD/out/multiqc/Reg1/after_merging " 2>&1 | tee \
	"$WD/log/multiqc/Reg1/after_merging/Report.log"

# Crear directorios
mkdir -p $WD/out/assembly/Reg1

# Indexar el genoma de referencia
bwa index $WD/res/Qbeta.fasta

# Alinear
find $WD/out/flash2/Reg1 -name "*.extendedFrags.fastq.gz" | sed \
	"s/.extendedFrags.fastq.gz$//" | parallel -j $(nproc) "bwa mem \
	$WD/res/Qbeta.fasta \
	$WD/out/flash2/Reg1/{/.}.extendedFrags.fastq.gz > \
	$WD/out/assembly/Reg1/{/.}.sam"

# Ordenar archivos SAM (human readable) y convertirlos en BAM (binarios)
find $WD/out/assembly/Reg1 -name "*.sam" | sed "s/.sam$//" | parallel -j \
	$(nproc) "samtools sort -o $WD/out/assembly/Reg1/{/.}.bam \
	$WD/out/assembly/Reg1/{/.}.sam"

# Indexar archivos BAM
find $WD/out/assembly/Reg1 -name "*.bam" | sed "s/.bam$//" | parallel -j \
	$(nproc) "samtools index $WD/out/assembly/Reg1/{/.}.bam"

# Crear directorios
mkdir -p $WD/out/assembly_filter/Reg1 $WD/out/assembly_not_filtered/Reg1

# Filtrar secuencias que no se alinearon
find $WD/out/assembly/Reg1 -name "*.bam" | sed "s/.bam$//" | parallel -j \
	$(nproc) "samtools view -F 0x04 -L $WD/data/Reg1/QBeta_R1.bed -U \
	$WD/out/assembly_not_filtered/Reg1/{/.}_not_filtered.bam -b \
	$WD/out/assembly/Reg1/{/.}.bam > \
	$WD/out/assembly_filter/Reg1/{/.}_filter.bam"

# Ordenar archivos BAM
find $WD/out/assembly_filter/Reg1 -name "*.bam" | sed "s/.bam$//" | parallel \
	-j $(nproc) "samtools sort -o $WD/out/assembly_filter/Reg1/{/.}_sort.bam \
	$WD/out/assembly_filter/Reg1/{/.}.bam"

# Indexar archivos BAM
find $WD/out/assembly_filter/Reg1 -name "*_sort.bam" | sed "s/_sort.bam$//" | parallel \
	-j $(nproc) "samtools index \
	$WD/out/assembly_filter/Reg1/{/.}_sort.bam"

# Crear directorios
mkdir -p $WD/out/Seqs_filter_length/Reg1

# Filtrar lecturas
for f in $WD/out/assembly_filter/Reg1/*_sort.bam; do samtools view -h \
	$f | awk '$1 ~ "^@" || length($10) == 271 && $6 !~ "I|D"' | samtools view \
	-b - > $WD/out/Seqs_filter_length/Reg1/$(basename $f .bam)_filter_length.bam; done;

# Crear directorios
mkdir -p $WD/out/seqs_indels/Reg1

# Filtrar lecturas con indels y gaurdar en un archivo
for f in $WD/out/assembly_filter/Reg1/*_sort.bam; do samtools view -h \
	$f | awk '$1 ~ "^@" || $6 ~ "I|D"' | samtools view -b - > \
	$WD/out/seqs_indels/Reg1/$(basename $f .bam)_indels.bam; done;

# Crear directorio
mkdir -p $WD/out/seqs_fasta/Reg1

# Convertir bam en fasta
find $WD/out/Seqs_filter_length/Reg1 -name "*.bam" | sed "s/.bam$//" | parallel \
	-j $(nproc) "samtools fasta \
	$WD/out/Seqs_filter_length/Reg1/{/.}.bam > \
	$WD/out/seqs_fasta/Reg1/{/.}.fasta"

# Crear directorios
mkdir -p $WD/out/seqs_fasta_collapsed/Reg1 $WD/log/seqs_fasta_collapsed/Reg1

# Colapsar lecturas
find $WD/out/seqs_fasta/Reg1 -name "*.fasta" | sed "s/.fasta$//" | parallel \
	-j $(nproc) "fastx_collapser -v -i \
	$WD/out/seqs_fasta/Reg1/{/.}.fasta -o \
	$WD/out/seqs_fasta_collapsed/Reg1/{/.}_collapsed.fasta 2>&1 | tee \
	$WD/log/seqs_fasta_collapsed/Reg1/{/.}.log"