#!/bin/bash
#Usage ShortReadAssembly.sh <input_directory_of fastq files>

#Check dependencies
echo  "Checking dependencies... "
for name in trimmomatic bwa samtools spades fastqc multiqc quast
do
  [[ $(which $name 2>/dev/null) ]] || { echo "$(tput setaf 1)$name$(tput sgr0) needs to be installed. Use 'conda install -c bioconda $name'";deps=1; }
done
[[ $deps -ne 1 ]] && echo "OK" || { echo -e "\n $(tput setaf 1)[ERROR] Install the above and rerun this script$(tput sgr0) \n You can find more information from github.com/tuc289/GABI";exit 1; }

#Illumina sequencing file uaually have two different system of file name
#If sequences were just out from the sequencer, File name contains "SampleName_S??_R[1,2]_001.fastq.gz", Next few lines will change the name as "SampleName_[1,2].fastq.gz"

cd $1

for f in *_R1_001.fastq.gz
do
	mv $f ${f%_S*_R1_001.fastq.gz}_1.fastq.gz
done

for f in *_R2_001.fastq.gz
do
	mv $f ${f%_S*_R2_001.fastq.gz}_2.fastq.gz
done

#If you download sequences from NCBi and dump fastq file (as fastq.gz) your filename should be already "SampleName_[1,2].fastq.gz"

#Running fastqc to check the reads quality
for file in *.fastq.gz
do
	if [ -d "${file%.fastq.gz}_fastqc" ]
	then
		echo "skip ${file}"
		continue
	fi
	echo "Running fastqc ${file}"
	fastqc $file -o ./
done

#echo "Running multiqc to combine fastqc reports"
#multiqc .

#Running trimmomatic to trim adapters and remove disqualified reads/bases

for f in *_1.fastq.gz
do
	if [ -f "${f%_1.fastq.gz}_1.trimmedP.fastq.gz}" ]
	then
	echo "skip ${f%_1.fastq.gz}"
	continue
	fi
	echo "Running Trimmomatic for ${f%_1.fastq.gz}"
	trimmomatic PE -threads 4 -phred33 $f ${f%_1.fastq.gz}_2.fastq.gz ${f%_1.fastq.gz}_1.trimmedP.fastq.gz ${f%_1.fastq.gz}_1.trimmedS.fastq.gz ${f%_1.fastq.gz}_2.trimmedP.fastq.gz ${f%_1.fastq.gz}_2.trimmedS.fastq.gz ILLUMINACLIP:/gpfs/group/jzk303/default/data/tuc289/rhAMR/amrplusplus_v2/data/adapters/nextera.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

#Running SPAdes to assemble genome

for f in *_1.trimmedP.fastq.gz
do
	if [ -d "${f%_1.trimmedP.fastq.gz}" ]
then
	echo 'skip ${f}'
	continue
fi
echo "Assemble ${f%_1.trimmedP.fastq.gz}"
spades.py -k 99,127 --isolate -1 $f -2 ${f%_1.trimmedP.fastq.gz}_2.trimmedP.fastq.gz -o ${f%_1.trimmedP.fastq.gz} -t 20 -m 64
done

mkdir contigs

for f in *_1.trimmedP.fastq.gz
do
	cd ${f%_1.trimmedP.fastq.gz}
	cat contigs.fasta > ${f%_1.trimmedP.fastq.gz}_contigs.fasta
	cp ${f%_1.trimmedP.fastq.gz}_contigs.fasta ../contigs
	cd ..
done

#Running quast to evaluate assembled genome
cd contigs
mkdir quast_results
quast *.fasta -o ../quast_results --min-contig 1000


#Calculate average coverage of the genomes
for f in *_contigs.fasta
do
	echo "Indexing $f"
	bwa index $f

	echo "Mapping reads to $f"
	bwa mem -t 20 $f ../${f%_contigs.fasta}_1.trimmedP.fastq.gz ../${f%_contigs.fasta}_2.trimmedP.fastq.gz -o ${f%_contigs.fasta}.sam
	echo "SAM created"

	echo "Converting SAM to BAM"
	samtools view -Sb ${f%_contigs.fasta}.sam -o ${f%_contigs.fasta}.bam
	echo "BAM created"
	rm *.sam

	echo "Sorting BAM file"
	samtools sort ${f%_contigs.fasta}.bam -o ${f%_contigs.fasta}_sorted.bam
	echo "Finished sorting"

	samtools index ${f%_contigs.fasta}_sorted.bam
	echo "Index complete"

	echo "Calcaulating average genome coverage"
	X=$(samtools depth ${f%_contigs.fasta}_sorted.bam | awk '{sum+=$3} END { print sum/NR}')
	echo "${f%_contigs.fasta}"
	echo "$X"
	echo "${f%_contigs.fasta} $X" > ${f}_average_coverage.txt
	rm *.bam
done

cat *_average_coverage.txt > average_coverage.txt

#File re-location
mv *.txt ../
mv *.fasta ../
cd ..
rm -r contigs

mv quast_results 6_quast
mkdir 1_fastq_raw
mkdir 2_fastqc_results
mkdir 3_fastq_trimmed
mkdir 4_spades
mkdir 5_contigs
mkdir 7_average_coverage

for f in *_1.fastq.gz
do
	mv ${f%_1.fastq.gz} 4_spades
done

mv *.trimmed* 3_fastq_trimmed
mv *.fastq.gz 1_fastq_raw
mv *_fastqc* 2_fastqc_results
mv *.fasta 5_contigs
mv *.txt 7_average_coverage

#Lastly, making combined report 
cp 6_quast/transposed_report.tsv ./quast_report.tsv
cp 7_average_coverage/average_coverage.txt ./


echo "Pipeline completed"