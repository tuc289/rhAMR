#Trimming reads based on the reads length and minimum length
for f in *_R1_001.fastq.gz
do
trimmomatic PE -threads 8 -phred33 -trimlog log $f ${f%_R1_001.fastq.gz}_R2_001.fastq.gz ${f%_R1_001.fastq.gz}_R1_001.trimmedP.fastq.gz ${f%_R1_001.fastq.gz}_R1_001.trimmedS.fastq.gz ${f%_R1_001.fastq.gz}_R2_001.trimmedP.fastq.gz ${f%_R1_001.fastq.gz}_R2_001.trimmedS.fastq.gz LEADING:3 TRAILING:3 MINLEN:50;
done
 
#merge paired-end reads in FLASH (output: merged_*.fastq)
for f in *.trimmedP.fastq.gz
do
flash $f ${f%_R1_001.trimmedP.fastq.gz}_R2_001.trimmedP.fastq.gz -M 150 --output-prefix merged_${f%_L001_R1_001.trimmedP.fastq.gz}
done
 
#convert fastq to fasta
for f in *.extendedFrags.fastq 
do
seqkit fq2fa $f --out-file ${f%.fastq}.fasta
done
 
#make custom blast database using Megares fasta file
makeblastdb -dbtype nucl -in full_sequence.fasta -out blastdb/mydb
blastdbcmd -db blastdb/mydb -info
 
#Local alignment with 100% identity
for f in merged_*.fasta
do
blastn -db blastdb/mydb -query $f -outfmt '6 qacc sacc pident length qcovs' -out blasted_${f%.fasta}.txt -num_threads 8 -perc_identity 100
done
