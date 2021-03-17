#!/bin/bash -x
#SBATCH -A SNIC2020-5-300
#SBATCH --output=/proj/nobackup/snic2019-35-62/gpozzati/job_out/gr-%A_%a.out
#SBATCH --error=/proj/nobackup/snic2019-35-62/gpozzati/job_err/gr-%A_%a.err
#SBATCH --array=1-223
#SBATCH -c 7
#SBATCH -t 01:00:00

list=$1
data=$2
out=$3
col=$4
offset=$5

pos=$(($SLURM_ARRAY_TASK_ID + $offset))
id=`tail -n+$pos $list | head -n 1`
echo 'Processing ' $id '...'

cd $data
ml singularity/3.5.3

if [[ ! -f $out ]]
then
    mkdir $out
fi


singularity exec --bind /proj/nobackup/snic2019-35-62/gpozzati/:/proj/nobackup/snic2019-35-62/gpozzati/ --nv /proj/nobackup/snic2019-35-62/gpozzati/sif/trfull.sif python /proj/nobackup/snic2019-35-62/gpozzati/BindingSite_docking_master/protocol_gramm.py -s1 processed_b4/$id'_u1.pdb' -s2 processed_b4/$id'_u2.pdb' -o $out/$id -g gramm-res/$id'_u1-'$id'_u2.res' -i1 formatted_labels/$id'_u1.site' -i2 formatted_labels/$id'_u2.site' -c $col

for n in {1..10}; do cat processed_b4/$id'_u1.pdb' | grep -v 'END' > $out/$id'_'$n'-tmp'; echo '' >> $out/$id'_'$n'-tmp';  cat $out/$id'_'$n'.pdb' >> $out/$id'_'$n'-tmp'; mv $out/$id'_'$n'-tmp' $out/$id'_'$n'.pdb'; done
