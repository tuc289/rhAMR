## Change "main_AmrPlusPlus_v2.nf" file
## Variable "java" should be deleted (line 74) "${java} -jar" should be deleted
 
## modify "config/local.config" file to indicate file path of dependencies
 
 
# AMRplusplus (megares)
export PATH=$PATH:/gpfs/group/jzk303/default/data/tuc289/rhAMR/amrplusplus_v2/bin:/storage/work/tuc289/miniconda3/envs/thesis/bin
 
/gpfs/group/jzk303/default/data/tuc289/rhAMR/amrplusplus_v2/nextflow run main_AmrPlusPlus_v2.nf -profile local --reads "/gpfs/group/jzk303/default/data/tuc289/rhAMR/illumina_1st/*_R{1,2}_001.fastq.gz"  --output output --threshold 80 
#Threshold can be changed (evantually, it needs to be set as "0", see below explanation from results)

rhAMR_illumina_1st

