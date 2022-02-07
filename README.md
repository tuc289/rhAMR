# rhAMR analysis archive 

This Repository will be used to archive essential script for rhAMR sequencing analysis

## Getting start

### Table of contents

#### 1. Whole genome sequencing assembly and AMR gene alignment to Megares (v2) database

Genome assemblies were done by running i) Trimmomatic 2) SPAdes and evaluated by i) fastqc ii) Quast iii) BWA, samtools (for calculating average coverage), followed by local pipeline ([ShortReadAssembly.sh](https://github.com/tuc289/rhAMR/ShortReadAssembly.sh)
). 

AMR gene alignmnet were done by Abricate (in local computer)

```
# download latest version of the database (Megares v2)
abricate-get_db --db megares

# Run abricate
abricate *_contigs.fasta > results.tab
abricate --summary results.tab > summary.tab
```
* in final summary.tab file, An absent gene is denoted . and a present gene is represented by its '%COVERAGE'.

#### 2. Human fecal metagenomics samples processing and AMR gene alignment (MG-Rast)
#### 3. rhAMR data analysis by AMRplusplus2
#### 4. Visualization of the results
