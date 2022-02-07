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

Raw sequences were uploadted to MG-Rast for AMR gene detection (This step can be done using AMRplusplus2, is it better to switch or not?)

#### 3. rhAMR data analysis by AMRplusplus2

##### Installation and configuration of AMRplusplus2 for running it in the Roar #####

```
#Install nextflow
$ wget -qO- https://get.nextflow.io | bash

#Clone pipeline source code
$ git clone https://github.com/meglab-metagenomics/amrplusplus_v2.git .

#Add current directory to $PATH
$ export PATH=$PATH:$(pwd)
```

1. In case you installed all the dependencies in local directory, you should add all of them into your PATH variable, or provide absolute path in the local configuration file **(config/local.config). 

2. If you are using **cond** , it could be easier to configure. First, check if all the dependencies were installed properly. Now, open the configure file **(config/local.config)** and edit as below.  

```
// The location of each dependency binary needs to be specified here.
// The examples listed below are assuming the tools are already in your $PATH, however,
// the absolute path to each tool can be entered individually.
env {
    /* These following tools are required to run AmrPlusPlus*/
    JAVA = "java" #DELETE THIS LINE
    TRIMMOMATIC = "trimmomatic-0.36.jar" #CHANGE IT AS "trimmomatic"
    PYTHON3 = "python3"
    BWA = "bwa"
    SAMTOOLS = "samtools"
    BEDTOOLS = 	"bedtools"
    RESISTOME = 	"resistome"
    RAREFACTION = 	"rarefaction"
    SNPFINDER = 	"snpfinder"
    FREEBAYES = "freebayes"
    /* These next tools are optional depending on which analyses you want to run */
    KRAKEN2 = "kraken2"
    RGI = "rgi"
    DIAMOND = "diamond"
}

process {
    cpus = 4                     // The maximum amount of CPUs to use
    disk = '125 GB'              // The maximum amount of disk space a single process is allowed to use
    //errorStrategy = 'ignore'     // Ignore process errors
    executor = 'local'           // The type of system the processes are being run on (do not modify this)
    maxForks = 1                 // The maximum number of forks a single process is allowed to spawn
    memory = '8 GB'              // The maximum amount of memory a single process is allowed to use
}
```

Now, open "main_AmrPlusPlus_v2.nf", and delete ${JAVA} from line 74.

##### Running AmrPlusPlus2

```
nextflow run main_AmrPlusPlus_v2.nf -profile local 
                                     --reads "<directory to raw reads>/*_R{1,2}_001.fastq.gz" 
                                     --output output --threshold 0
```
##### Check the results



#### 4. Visualization of the results
