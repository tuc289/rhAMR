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

2. If you are using **conda** , it could be easier to configure. First, check if all the dependencies were installed properly. Now, open the configure file **(config/local.config)** and edit as below.  

```
// The location of each dependency binary needs to be specified here.
// The examples listed below are assuming the tools are already in your $PATH, however,
// the absolute path to each tool can be entered individually.
env {
    /* These following tools are required to run AmrPlusPlus*/
    JAVA = "java" ###############################DELETE THIS LINE###########################
    TRIMMOMATIC = "trimmomatic-0.36.jar" ######################CHANGE IT AS "trimmomatic"##################
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
The program will generate multiple directories after it is completed
```
Some notable results 
RunResistome - .type.tsv, .class.tsv, .gene.tsv, .mechanism.tsv for individual samples (#of Hits)
ResistomeResults - Matrix output from resistome (AMR genes) with #of Hits and database match (AMR_analytic_matrix.csv)
AlignToAMR - .bam and .sam file abouot alignment results from reads to the database
```

#### 4. Visualization of the results (in progress)

First, AMR_analytic_matrix.csv should be modified for loading it to R (i.e., remove unwanted column, seperate "|" deliminated column into multiple for better classification

The, save it as .csv file and load it into R, and manipulate it for plotting
```
for_R <- read.table("~/Downloads/REsistomeResults/AMR_analytic_matrix.csv", sep=",", header=T, row.names=1)
str(for_R)
colnames(for_R) <- c("Type", "Class", "Mechanism", "Gene", "p0_1:100", "p0_1:1000", 
                     "p0_1:10000", "p2_1:1", "p2_1:100", "p2_1:1000", "p2_1:10000")
for_plot <- for_R[,4:11]
```

Now, extract targeted genes for LOD analysis
```
p0_gene_list <- c("ACRA", "ACRB", "ACRD", "ACRE", "ACRF", "ACRS", "AMPH", "ASMA", "BACA", 
                  "BAER", "BAES", "BLAEC", "CATA", "CPXAR", "CRP", "DFRA", "EMRA", "EMRB", "EMRD",
                  "EMRK", "EMRR", "EMRY", "EPTA", "EVGS", "GADW", "HNS", "KDPE", "KPNO", "MARA", "MARR", 
                  "MDFA", "MDTA", "MDTB", "MDTC", "MDTE", "MDTF", "MDTG", "MDTH", "MDTI", "MDTJ",
                  "MDTK", "MDTN", "MDTO", "MDTP", "MPHA", "MPHB", "MSBA", "MVRC", "PBP4B", "PMRF", 
                  "ROBA", "SOXS", "TETB")
p2_gene_list <- c("ACRB", "DFRA", "GADW", "MDTK")

p0_plot <- filter(for_plot, Gene %in% p0_gene_list)
p2_plot <- filter(for_plot, Gene %in% p2_gene_list)
p2_plot_genes <- p2_plot[1]

p0_plot <- p0_plot[1:4]
p2_plot <- p2_plot[5:8]
p2_plot <- cbind(p2_plot_genes, p2_plot)
```

Melt the dataframe into multiple rows (this will help plotting bar_plot in each combination of sample X gene)
```
library(reshape2)
melted_p0 <- melt(p0_plot)
melted_p2 <- melt(p2_plot)

head(melted_p0)
melted_p0
head(melted_p2)
melted_p2
```

Creating multiple bar plots and save it to the computer
```
p0_plot = ggplot(melted_p0, aes(x = variable, fill = Gene, y = value)) + 
  #facet_grid(vars(Gene), scales = "free") +
  geom_bar(stat = "identity") + 
  facet_wrap(~ Gene, nrow=10, ncol=10, scales = "free") + 
  #theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1, face= "bold"), 
  #      axis.title.y = element_text(size = 16, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
  #      legend.text = element_text(size = 12, face = "bold", colour = "black"), 
  #      axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
  scale_y_continuous(trans='log10') +
  guides(fill="none") +
  xlab("Ratio")+
  ylab("Count")
p0_plot

ggsave("p0_plots.tiff", plot=p0_plot, device=tiff, width=30, height=15, units=c("in"))
```
Do the same thing with different object name for other panel and/or samples

***To-do list***

Running regression analysis on each gene - list all R2 and p values

re-creating line plot with multiple lines (indicating each gene)

set-up the cut-off for detecting right LOD




