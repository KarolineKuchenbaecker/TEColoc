#!/bin/bash

TraitName=$1
# should also be the name of the phenotype column for study 2 (raw data)
chr=$2
rs=$3
bp=$4
range=$5
################

# select population for study1 (the sample set with summary statistics)
pop=EUR

###################
# set paths correctly
# folder for all the temporary files (will be created if it does not exist)
tmp=~tmp_$TraitName
# file with the summary stats for study 1 
s1sumstats=~/res_${TraitName}.txt.gz
# folder where extracted summary stats should be stored
s1tmpdir=~/sumstats/
# where the results files should be written to
results=~/${TraitName}.${chr}.${bp}.jlim.out
# plink file (base name) with raw genotype data for study 2. assuming separate files for each chr ending on the chr (which is added automatically)
s2plink=~/myplink
# phenotype data for study 2
s2phenofile=~/myplink_phenos
# folder where jlim is located 
jlimf=~/bin/
# folder with the 1000 genomes v3 data as vcfs
tgf=~/tg/


#########################
### for basic run no need to modify from here on 
let start=$bp-$range*1000
let end=$bp+$range*1000
cfg=${tmp}/locus.${chr}.${start}.${end}/jlim.cfg.${TraitName}.tsv

if [ ! -d ${tmp} ]; then
  mkdir -p ${tmp};
fi
mkdir  ${tmp}/locus.${chr}.${start}.${end}

# extract summary stats
echo "SNP CHR BP P" > ${s1tmpdir}${TraitName}.${chr}.${start}.${end}.txt
tabix $s1sumstats ${chr}:${start}-${end} | awk -v OFS=' ' '{print $1,$3,$5,$16}' >> ${s1tmpdir}${TraitName}.${chr}.${start}.${end}.txt
grep $rs ${s1tmpdir}${TraitName}.${chr}.${start}.${end}.txt


# prep ukhls geno 
plink --bfile ${s2plink}${chr} --recode --geno 0.0 --chr $chr --from-bp $start --to-bp $end --out ${tmp}/locus.${chr}.${start}.${end}/$TraitName --maf 0.005
plink --file ${tmp}/locus.${chr}.${start}.${end}/$TraitName  --pheno ${s2phenofile} --pheno-name ${TraitName}  --mperm 1000 --linear --mperm-save-all --out  ${tmp}/locus.${chr}.${start}.${end}/${TraitName}.nix 

# # # # implement a step where it checks the snp overlap between studies.
# only proceed of the overlap is large enough 
tr1extract=${s1tmpdir}${TraitName}.${chr}.${start}.${end}.txt
tr2=${tmp}/locus.${chr}.${start}.${end}/${TraitName}.map
#comm -12 <(cut -f1 $tr1 -d' ' | tail -n +2 | sort ) <(cut -f2 $tr2 | sort ) 
overlap="$(comm -12 <(cut -f1 $tr1extract -d' ' | tail -n +2 | sort ) <(cut -f2 $tr2 | sort ) | wc -w)"
n1="$(cut -f1 $tr1extract -d' ' | tail -n +2 | wc -l )"
n2="$(cut -f2 $tr2 | wc -l )"
echo $overlap $n1 $n2

if [ $overlap -ge 10 ]
then

gzip ${tmp}/locus.${chr}.${start}.${end}/*

echo -e "CHR\tSNP\tBP\tSTARTBP\tENDBP" > ${tmp}/locus.${chr}.${start}.${end}/indexSNP.tsv
echo -e "${chr}\t${rs}\t${bp}\t${start}\t$end" >> ${tmp}/locus.${chr}.${start}.${end}/indexSNP.tsv

${jlimf}jlim/bin/fetch.refld0.${pop}.pl ${tgf} ${tmp}/locus.${chr}.${start}.${end}/indexSNP.tsv ${tmp}/locus.${chr}.${start}.${end}/

${jlimf}jlim/bin/jlim_gencfg.sh --tr1-name ${TraitName} --tr1-dir $s1tmpdir --tr2-dir ${tmp} --idxSNP-file ${tmp}/locus.${chr}.${start}.${end}/indexSNP.tsv --refld-dir ${tmp}/locus.${chr}.${start}.${end}/ --out $cfg --p-tr2-cutoff 0.01

${jlimf}jlim/bin/run_jlim.sh $cfg 0.8 $results

if [ -f $results ]
then
cut -f1-7,12 $cfg > ${rs}tmp1
echo -e  "overlap\tn1\tn2" > ${rs}tmp2
echo -e "${overlap}\t$n1\t$n2" >> ${rs}tmp2
paste -d'\t' ${rs}tmp1 ${rs}tmp2 > ${rs}tmp3 
paste -d'\t' ${rs}tmp3 $results > ${results}.full 
rm ${rs}tmp1 ${rs}tmp2 ${rs}tmp3 
cat ${results}.full 
else
	echo "analysis did not produce results"
fi

else
	echo -e "${TraitName}\t${chr}\t${rs}\t${bp}\t${start}\t${end}\t${overlap}\t${n1}\t${n2}" >> ${tmp}/low_overlap.txt
fi

######################
