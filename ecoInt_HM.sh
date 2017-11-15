#!/bin/bash
#SBATCH --job-name=AGecoint
#SBATCH --output=/work_zfs2/sukem091/AGORA-interactions-MTF2/EcoInt.log
#SBATCH --error=/work_zfs2/sukem091/AGORA-interactions-MTF2/EcoInt.err
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=90G
#SBATCH --time=36:00:00
#SBATCH --partition=msb
#SBATCH --qos=msb

scriptpath="/work_zfs2/sukem091/AGORA-interactions-MTF2"

models="agora_HumanMilkDietMM"
swcorr="SWcorr"
minGrowth="nminGrowth"
cores="30"


cd $TMPDIR
echo $TMPDIR
date

cp $scriptpath/EcoInt_cluster.R .
cp $scriptpath/join_mult_models.R .
cp $scriptpath/coupling.R .
cp $scriptpath/addMultiReact.R .
cp $scriptpath/get_metabolic_interchange.R .
cp /work_zfs2/sukem091/Resources/AGORA/correct_common_errors_agora2.R .

cp /work_zfs2/sukem091/Resources/AGORA/$models.RDS .
mkdir results

module load R-3.2.0
Rscript EcoInt_cluster.R $models $swcorr $minGrowth $cores

cp results/* $scriptpath/results/
