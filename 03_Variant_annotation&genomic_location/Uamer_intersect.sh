#!/bin/bash
#SBATCH -N 1 # number of nodes
#SBATCH -n 16 # number of cores
#SBATCH --mem=48G
#SBATCH --job-name="cactus_halLiftover"
#SBATCH --mail-type=END
#SBATCH --partition=computeq
#SBATCH --mail-user=hrclndnn@memphis.edu
#SBATCH -a 0-18

module load bcftools
module load bedtools
module load cactus
module load htslib


#SnpEff annotated vcf: /home/hrclndnn/snpEff/Uamer/Uamer.ann.vcf.gz
#formatted ROH bed files: /home/hrclndnn/snpEff/Uamer/liftover/${file}/${file}_ROH.bed


files=(AK17023 
AK17047 
AK17117 
AZ12 
AK18242 
ID10 
LA366 
LA371 
LAT594 
MI334 
MI335 
MN6083 
MS3783 
NC00417 
NC056 
NVb84 
NVb99 
NVg5 
WV17)

file=${files[$SLURM_ARRAY_TASK_ID]}

#00 set up directories
cd ./maf55
mkdir ./${file}

#01 Create vcfs for each individual
bcftools view -Oz -o ./${file}/${file}.maf55.vcf.gz -s ${file} -q 0.55 /home/hrclndnn/snpEff/Uamer/liftover/renamed_Uamer.ann.vcf.gz 

#02 Create filtered bedfiles for each individual from annotations
#pull out scaffold name, start position, create end position, and include annotation field as unique ID for variants of interest 
bcftools view -H ./${file}/${file}.maf55.vcf.gz | grep frameshift_variant | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$8}' > ./${file}/${file}.maf55_frameshift_variant.out
bcftools view -H ./${file}/${file}.maf55.vcf.gz | grep missense_variant | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$8}' >./${file}/${file}.maf55_missense_variant.out
bcftools view -H ./${file}/${file}.maf55.vcf.gz | grep initiator_codon_variant | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$8}' >./${file}/${file}.maf55_initiator_codon_variant.out
bcftools view -H ./${file}/${file}.maf55.vcf.gz | grep stop_retained_variant | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$8}' >./${file}/${file}.maf55_stop_retained_variant.out
bcftools view -H ./${file}/${file}.maf55.vcf.gz | grep rare_amino_acid_variant | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$8}' >./${file}/${file}.maf55_rare_amino_acid_variant.out
bcftools view -H ./${file}/${file}.maf55.vcf.gz | grep splice_acceptor_variant | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$8}' >./${file}/${file}.maf55_splice_acceptor_variant.out
bcftools view -H ./${file}/${file}.maf55.vcf.gz | grep splice_donor_variant | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$8}' >./${file}/${file}.maf55_splice_donor_variant.out
bcftools view -H ./${file}/${file}.maf55.vcf.gz | grep stop_lost | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$8}' >./${file}/${file}.maf55_stop_lost.out
bcftools view -H ./${file}/${file}.maf55.vcf.gz | grep 5_prime_UTR_premature | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$8}' >./${file}/${file}.maf55_5_prime_UTR_premature.out
bcftools view -H ./${file}/${file}.maf55.vcf.gz | grep start_lost | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$8}' >./${file}/${file}.maf55_start_lost.out
bcftools view -H ./${file}/${file}.maf55.vcf.gz | grep stop_gained | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$8}' >./${file}/${file}.maf55_stop_gained.out
bcftools view -H ./${file}/${file}.maf55.vcf.gz | grep synonymous_variant | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$8}' >./${file}/${file}.maf55_synonymous_variant.out
bcftools view -H ./${file}/${file}.maf55.vcf.gz | grep start_retained | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$8}' >./${file}/${file}.maf55_start_retained.out
bcftools view -H ./${file}/${file}.maf55.vcf.gz | grep stop_retained_variant | awk 'BEGIN {OFS="\t"} {print $1,$2,$2+1,$8}' >./${file}/${file}.maf55_stop_retained_variant.out

#03 Concatenate individuals' variants into pseudo-bed file format
#maf filtered
cat ./${file}/${file}.maf55_frameshift_variant.out ./${file}/${file}.maf55_missense_variant.out ./${file}/${file}.maf55_initiator_codon_variant.out ./${file}/${file}.maf55_stop_retained_variant.out ./${file}/${file}.maf55_rare_amino_acid_variant.out ./${file}/${file}.maf55_splice_acceptor_variant.out ./${file}/${file}.maf55_splice_donor_variant.out ./${file}/${file}.maf55_stop_lost.out ./${file}/${file}.maf55_5_prime_UTR_premature.out ./${file}/${file}.maf55_start_lost.out ./${file}/${file}.maf55_stop_gained.out ./${file}/${file}.maf55_synonymous_variant.out ./${file}/${file}.maf55_start_retained.out ./${file}/${file}.maf55_stop_retained_variant.out > ./${file}/${file}.maf55_all_variants.out

#04 Remove duplicate lines & sort
sort -u ./${file}/${file}.maf55_all_variants.out > ./${file}/${file}.maf55_allVar_dupRemoved.txt
#count variants

#05 Use cactus to convert positions from NCBI assembly to DNAzoo assembly
#module load cactus
halLiftover /home/mdpllard/Cactus/heather_bears/steps-output/bears_alignment.hal Ursus_americanus_NCBI /home/hrclndnn/snpEff/Uamer/maf55/${file}/${file}.maf55_allVar_dupRemoved.txt Ursus_americanus_DNAZoo /home/hrclndnn/snpEff/Uamer/maf55/${file}/${file}.maf55_liftover.bed
#remove duplicates
cat ./${file}/${file}.maf55_liftover.bed | sort -k4,4 | uniq -u -f3 > ./${file}/${file}.maf55_liftover_noDups.bed

#06 remove all but scaffolds 2-37 for maf filtered:
ls -l | grep 'HiC_scaffold_2[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold2.bed
ls -l | grep 'HiC_scaffold_3[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold3.bed
ls -l | grep 'HiC_scaffold_4[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold4.bed
ls -l | grep 'HiC_scaffold_5[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold5.bed
ls -l | grep 'HiC_scaffold_6[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold6.bed
ls -l | grep 'HiC_scaffold_7[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold7.bed
ls -l | grep 'HiC_scaffold_8[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold8.bed
ls -l | grep 'HiC_scaffold_9[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold9.bed
ls -l | grep 'HiC_scaffold_10[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold10.bed
ls -l | grep 'HiC_scaffold_11[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold11.bed
ls -l | grep 'HiC_scaffold_12[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold12.bed
ls -l | grep 'HiC_scaffold_13[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold13.bed
ls -l | grep 'HiC_scaffold_14[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold14.bed
ls -l | grep 'HiC_scaffold_15[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold15.bed
ls -l | grep 'HiC_scaffold_16[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold16.bed
ls -l | grep 'HiC_scaffold_17[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold17.bed
ls -l | grep 'HiC_scaffold_18[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold18.bed
ls -l | grep 'HiC_scaffold_19[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold19.bed
ls -l | grep 'HiC_scaffold_20[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold20.bed
ls -l | grep 'HiC_scaffold_21[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold21.bed
ls -l | grep 'HiC_scaffold_22[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold22.bed
ls -l | grep 'HiC_scaffold_23[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold23.bed
ls -l | grep 'HiC_scaffold_24[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold24.bed
ls -l | grep 'HiC_scaffold_25[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold25.bed
ls -l | grep 'HiC_scaffold_26[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold26.bed
ls -l | grep 'HiC_scaffold_27[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold27.bed
ls -l | grep 'HiC_scaffold_28[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold28.bed
ls -l | grep 'HiC_scaffold_29[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold29.bed
ls -l | grep 'HiC_scaffold_30[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold30.bed
ls -l | grep 'HiC_scaffold_31[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold31.bed
ls -l | grep 'HiC_scaffold_32[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold32.bed
ls -l | grep 'HiC_scaffold_33[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold33.bed
ls -l | grep 'HiC_scaffold_34[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold34.bed
ls -l | grep 'HiC_scaffold_35[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold35.bed
ls -l | grep 'HiC_scaffold_36[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold36.bed
ls -l | grep 'HiC_scaffold_37[[:blank:]]' ./${file}/${file}.maf55_liftover_noDups.bed > ./${file}/${file}.maf55_scaffold37.bed

#concatenate scaffolds into one bed file
cat ./${file}/${file}.maf55_scaffold2.bed ./${file}/${file}.maf55_scaffold3.bed ./${file}/${file}.maf55_scaffold4.bed ./${file}/${file}.maf55_scaffold5.bed ./${file}/${file}.maf55_scaffold6.bed ./${file}/${file}.maf55_scaffold7.bed ./${file}/${file}.maf55_scaffold8.bed ./${file}/${file}.maf55_scaffold9.bed ./${file}/${file}.maf55_scaffold10.bed ./${file}/${file}.maf55_scaffold11.bed ./${file}/${file}.maf55_scaffold12.bed ./${file}/${file}.maf55_scaffold13.bed ./${file}/${file}.maf55_scaffold14.bed ./${file}/${file}.maf55_scaffold15.bed ./${file}/${file}.maf55_scaffold16.bed ./${file}/${file}.maf55_scaffold17.bed ./${file}/${file}.maf55_scaffold18.bed ./${file}/${file}.maf55_scaffold19.bed ./${file}/${file}.maf55_scaffold20.bed ./${file}/${file}.maf55_scaffold21.bed ./${file}/${file}.maf55_scaffold22.bed ./${file}/${file}.maf55_scaffold23.bed ./${file}/${file}.maf55_scaffold24.bed ./${file}/${file}.maf55_scaffold25.bed ./${file}/${file}.maf55_scaffold26.bed ./${file}/${file}.maf55_scaffold27.bed ./${file}/${file}.maf55_scaffold28.bed ./${file}/${file}.maf55_scaffold29.bed ./${file}/${file}.maf55_scaffold30.bed ./${file}/${file}.maf55_scaffold31.bed ./${file}/${file}.maf55_scaffold32.bed ./${file}/${file}.maf55_scaffold33.bed ./${file}/${file}.maf55_scaffold34.bed ./${file}/${file}.maf55_scaffold35.bed ./${file}/${file}.maf55_scaffold36.bed ./${file}/${file}.maf55_scaffold37.bed > ./${file}/${file}.maf55_scaffolds2-37.bed

#07 find intersection of bed files
#output both where intersection occurs and variant annotation info
#maf filtered
bedtools intersect -wa -wb -a /home/hrclndnn/snpEff/Uamer/liftover/${file}/${file}_ROH.bed -b ./${file}/${file}.maf55_scaffolds2-37.bed > ./${file}/${file}.maf55_intersect_rohVar.txt
