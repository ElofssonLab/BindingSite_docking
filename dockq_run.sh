#!/bin/bash -x
#SBATCH -A SNIC2020-5-300 
#SBATCH --output=/proj/nobackup/snic2019-35-62/gpozzati/job_out/dq-%A_%a.out
#SBATCH --error=/proj/nobackup/snic2019-35-62/gpozzati/job_err/dq-%A_%a.error
#SBATCH --array=1-223
#SBATCH -c 1
#SBATCH -t 01:00:00

list=$1
offset=$2
folder=$3

pos=$(($SLURM_ARRAY_TASK_ID + $offset))
id=`tail -n+$pos $list | head -n 1`
echo $id

ml singularity/3.5.3

for n in {1..10}; do
singularity exec --bind /proj/nobackup/snic2019-35-62/gpozzati/:/proj/nobackup/snic2019-35-62/gpozzati/ --nv /proj/nobackup/snic2019-35-62/gpozzati/sif/tf.sif /proj/nobackup/snic2019-35-62/gpozzati/programs/DockQ/DockQ.py -short /proj/nobackup/snic2019-35-62/gpozzati/BindingSite_docking_master/data/$folder/$id'_'${n}'.pdb' /proj/nobackup/snic2019-35-62/gpozzati/BindingSite_docking_master/data/processed_b4/$id'_b.pdb' | grep 'DockQ' >> /proj/nobackup/snic2019-35-62/gpozzati/BindingSite_docking_master/data/$folder/$id'_results'; done

