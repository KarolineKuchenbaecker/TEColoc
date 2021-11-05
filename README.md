# TAColoc
trans-ancestry colocalization 

Trans-ancestry colocalization (TAColoc) tests whether a specific locus has the same causal variant in two groups with different ancestry. It has been described in https://doi.org/10.1101/525170. It uses the Joint Likelihood Mapping (JLIM) model. JLIM was originally developed for colocalization testing of GWAS and eQTL signals. See http://genetics.bwh.harvard.edu/wiki/sunyaevlab/jlim or https://github.com/cotsapaslab/jlim for more details about JLIM.


# Requirements
It is necessary to install JLIM. For the first sample set, summary statistics are used. For that group LD is approximated by estimating it from a reference panel. In the current script, the 1000 Genomes Project data are used. Samples that match the study ancestry group for the summary statistics data set need to be selected from the 1000G set. We have provided a script for that. For the second ancestry group raw genotype data need to be available.
Additionally, plink and tabix are required.

Note that this code creates a certain folder structure. This originates from JLIM and cannot be modified without changing the JLIM code.


# Install JLIM 
git clone https://github.com/cotsapaslab/jlim.git
Rscript -e 'install.packages("getopt", repos="http://cran.r-project.org")' 
cd jlim
R-3 CMD INSTALL jlimR_1.0.2.tar.gz


# Select population with matching ancestry from 1000 Genomes Project
Use this script to modify the 1000 Genomes population that is used to calculate the LD scores. It should match the ancestry of the samples with summary statistics.
select_population_ldscores.sh

# Prepare input files
The raw data for the first ancestry group needs to be in plink binary format.
The file with the summary statistics is assumed to be a zipped (.gz) file that contains information for a SNP identifier, chromosome, position and p-value. It also needs to be indexed using tabix.

Index summary statistics file
resfile=/path/to/summarystats/file.txt 
bgzip $resfile 
tabix -c "SNP" -s 2 -b 3 -e 3 $resfile 
to test extract a given region
tabix ${resfile}.gz 16:56993320-56993326


# Run TAColoc
Directory paths for the raw data, summary stats and working directory need to be modified in the script manually. It can then be run for a specific locus via 
./TEColoc_v1.sh [PhenotypeName] [chromosome] [chromosome] [lead SNP] [position]
For example
./TEColoc_v1.sh LDL 19 rs4420638 45422946
