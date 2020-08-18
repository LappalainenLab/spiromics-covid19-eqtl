#!/bin/bash
#SBATCH --job-name=process                           # Job name
#SBATCH --mem=75gb                                  # Job memory request
##SBATCH --time=8:00:00                             # Time limit 4 hours
#SBATCH --output=log/process_chr_%A_%a.out           # By default stderr and stdout is redirected to the same file, use %j if not using task array
#SBATCH --error=log/process_chr_%A_%a.error           # Redirect stderr
#SBATCH --array=1-23                                 # Array range; 1-k, where k is the length of COV vector

## ---- Process TOPMed Freeze 9 unphased genotype data -------------------------
# Select only variants that PASS all the filters & select only 2,710 SPIROMICS samples
# Select common variants with MAF >= 1%
# Define variant IDs - chrXX_pos_ref_alt_b38

module load bcftools/1.9

# Path to data
dir='/path/to/freeze.9/minDP10'

# Define arguments
CHR=({1..22} X)

# Setup slurm task array
IDX=$(($SLURM_ARRAY_TASK_ID - 1))

# Note: Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF-to-BCF conversion.

bcftools view -Ou -f PASS --samples-file linking_files/master_with_resends_and_topmedids.cleaned_formatted.NWDID_present_in_freeze9_2710indiv.txt ${dir}/freeze.9.chr${CHR[${IDX}]}.pass_and_fail.gtonly.minDP10.bcf | bcftools view -Ou -i 'MAF>=0.01' | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT\_b38' -Ob -o genotype/tmp_chr${CHR[${IDX}]}.freeze9.pass_only.spiromics_2710samples.maf01.bcf
