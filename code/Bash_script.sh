### Download data
### METH dependent vs naive (PMID: 31605699, ENA: PRJEB33332)
cd /mnt/d/Projects/OxyMeth
mkdir data
cd data
mkdir PRJEB33332
cd PRJEB33332

#### Put fatsqs in PRJEB33332 folder, then QC
mkdir fastqc
for fname in *.fastq.gz
do
root=`basename $fname .fastq.gz`
echo "Doing $root"
fastqc ${root}.fastq.gz -o fastqc/ -t 12
done

# adapter trimming
mkdir trimmed
for fname in *_R1.fastq.gz
do
root=`basename $fname _R1.fastq.gz`
echo "Doing $root"
fastp -i ${root}_R1.fastq.gz -I ${root}_R2.fastq.gz -o trimmed/${root}_R1_trimmed.fastq.gz -O trimmed/${root}_R2_trimmed.fastq.gz
rm ${root}_R1.fastq.gz ${root}_R2.fastq.gz
done

### Create custom hybrid genome
# Rnor6

# Create Rnor 6 genome
# cd /mnt/f/genomes/Rat/FASTA
# cat /mnt/f/genomes/Rat/FASTA/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa | grep '>' | cut -f 1 -d ' ' | sed 's/>//g' > rn6_chromosomes.txt

# # now create custom fasta deleting unnecessary chromosome scaffolds
# seqtk subseq -l 60 Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa rn6_chromosomes.txt > rn6.fa 
# Create the hybrid genome
cd /mnt/f/genomes/RSEM/Rat/DNA/intermediate/
cat Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa AF324493.2.fasta > Rnor6_hybrid.fa
mv Rnor6_hybrid.fa > /mnt/f/genomes/Rat/Rnor6_HIV/Rnor6_hybrid.fa

# make hisat2 index
cd trimmed/
mkdir hisat_bams
# Hisat index
#index=/mnt/g/Sanna/Rat_Genome/hisat_index/Rat
# new hybrid index
index=/mnt/d/genomes/Rat/Rnor6_HIV/hisat/Rnor6

# Align with hisat2
for fname in *_R1_trimmed.fastq.gz
do
cd /mnt/d/Projects/OxyMeth/data/PRJEB33332/trimmed/
root=`basename $fname _R1_trimmed.fastq.gz`
echo "Doing $root"
hisat2 -x $index \
-p 12 \
-t \
-1 ${root}_R1_trimmed.fastq.gz \
-2 ${root}_R2_trimmed.fastq.gz \
| samtools view -bS - > hisat_bams/${root}.bam
rm ${root}_R1_trimmed.fastq.gz ${root}_R2_trimmed.fastq.gz
# Sort and index BAMs
cd /mnt/d/Projects/OxyMeth/data/PRJEB33332/trimmed/hisat_bams
samtools sort -m 2G -@ 3 -O BAM -o ${root}.sorted.bam ${root}.bam 
samtools index ${root}.sorted.bam
rm ${root}.bam
done

### Convert to counts with the subread package
cd /mnt/d/Projects/OxyMeth/data/PRJEB33332/trimmed/hisat_bams
mkdir rawcounts
#gtf=/mnt/g/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf
#featureCounts -T 7 -p -s 2 -t gene -g gene_id -a $gtf -o rawcounts/featureCounts.txt *sorted.bam

# new counts with hybrid genome
gtf=/mnt/d/genomes/Rat/Rnor6_HIV/Rnor6_hybrid_fixed.gtf
featureCounts -T 10 -p -s 2 -t gene -g gene_name -a $gtf -o rawcounts/featureCounts.txt *sorted.bam

# Unique mapping rates: 68.4,69.4,69.3,70.0,69.0,69.0,68.9,69.1,69.4,70.3,69.1,69.5,68.7,69.9,70.3
# mean 69.4

#### Second dataset, P35458399
cd /mnt/d/Projects/OxyMeth/data/P35458399

#### Put fatsqs in P35458399 folder, then
mkdir fastqc
for fname in *.fastq.gz
do
root=`basename $fname .fastq.gz`
echo "Doing $root"
fastqc ${root}.fastq.gz -o fastqc/ -t 12
done

# adapter trimming
mkdir trimmed
for fname in *_R1_001.fastq.gz
do
root=`basename $fname _R1_001.fastq.gz`
echo "Doing $root"
fastp -i ${root}_R1_001.fastq.gz -I ${root}_R2_001.fastq.gz -o trimmed/${root}_R1_trimmed.fastq.gz -O trimmed/${root}_R2_trimmed.fastq.gz
rm ${root}_R1_001.fastq.gz ${root}_R2_001.fastq.gz
done

mkdir hisat_bams
# Hisat index
#index=/mnt/g/Sanna/Rat_Genome/hisat_index/Rat
# new hybrid index
index=/mnt/d/genomes/Rat/Rnor6_HIV/hisat/Rnor6

# Align with hisat2
for fname in *_R1_trimmed.fastq.gz
do
cd /mnt/d/Projects/OxyMeth/data/P35458399/trimmed/
root=`basename $fname _R1_trimmed.fastq.gz`
echo "Doing $root"
hisat2 -x $index \
-p 12 \
-t \
-1 ${root}_R1_trimmed.fastq.gz \
-2 ${root}_R2_trimmed.fastq.gz \
| samtools view -bS - > hisat_bams/${root}.bam
rm ${root}_R1_trimmed.fastq.gz ${root}_R2_trimmed.fastq.gz
# Sort and index BAMs
cd /mnt/d/Projects/OxyMeth/data/P35458399/trimmed/hisat_bams
samtools sort -m 2G -@ 3 -O BAM -o ${root}.sorted.bam ${root}.bam 
samtools index ${root}.sorted.bam
rm ${root}.bam
done

### Convert to counts with the subread package
cd /mnt/d/Projects/OxyMeth/data/P35458399/trimmed/hisat_bams
mkdir rawcounts
#gtf=/mnt/g/Sanna/Rat_Genome/GTF/Rattus_norvegicus.Rnor_6.0.100.gtf
#featureCounts -T 7 -p -s 2 -t gene -g gene_id -a $gtf -o rawcounts/featureCounts.txt *sorted.bam

# new counts with hybrid genome
gtf=/mnt/d/genomes/Rat/Rnor6_HIV/Rnor6_hybrid_fixed.gtf
featureCounts -T 10 -p -s 2 -t gene -g gene_name -a $gtf -o rawcounts/featureCounts.txt *sorted.bam
