#!/bin/bash
#SBATCH -N 1 # number of nodes
#SBATCH -n 16 # number of cores
#SBATCH --mem=48G
#SBATCH --job-name="SnpEff"
#SBATCH --mail-type=END
#SBATCH --partition=computeq
#SBATCH --mail-user=hrclndnn@memphis.edu

module load snpeff

snpeff -v -Xmx4g "ASM334442v1.99" /home/hrclndnn/WGSUamer/09filterSNPs/Uamer.snp.filt.vcf > ./Uamer.ann.vcf
snpeff -v UrsMar_1.0.99 /home/hrclndnn/WGSUarct/Uarct_renamed238scaffolds.vcf > Uarctos.ann.vcf
snpeff -v UrsMar_1.0.99 /home/hrclndnn/WGSUmari/09filterSNPs/renamed238scaffolds.vcf > Umar.ann.vcf
