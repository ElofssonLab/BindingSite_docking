list=$1
out=$2
col=$3
cores=$4
folder=`pwd`

if [[ ! -f $out ]]
then
    mkdir results/$out
fi

while read id; do
    echo 'Processing ' $id '...'
    singularity exec --bind $folder:$folder --nv $folder/BSenv.sif python $folder/src/protocol_gramm.py -s1 $folder/data/processed_b4/$id'_u1.pdb' -s2 $folder/data/processed_b4/$id'_u2.pdb' -o $folder/results/$out/$id -g $folder/data/gramm-res/$id'_u1-'$id'_u2.res' -i1 $folder/data/formatted_labels/$id'_u1.site' -i2 $folder/data/formatted_labels/$id'_u2.site' -c $col -n $cores;
    for n in {1..10}; do 
	cat $folder/results/$out/$id'_r.pdb' | grep -v 'END' > $folder/results/$out/$id'_'$n'-tmp'; 
	cat $folder/results/$out/$id'_'$n'.pdb' >> $folder/results/$out/$id'_'$n'-tmp'; 
	mv $folder/results/$out/$id'_'$n'-tmp' $folder/results/$out/$id'_'$n'.pdb'; 
	singularity exec --bind $folder:$folder --nv $folder/BSenv.sif /DockQ/DockQ.py -short $folder/results/$out/$id'_'$n'.pdb' $folder/data/processed_b4/$id'_uc.pdb' | grep 'DockQ' >> $folder/results/$out/$id'_resultsu';
        singularity exec --bind $folder:$folder --nv $folder/BSenv.sif /DockQ/DockQ.py -short $folder/results/$out/$id'_'$n'.pdb' $folder/data/processed_b4/$id'_bc.pdb' | grep 'DockQ' >> $folder/results/$out/$id'_resultsb';
    done;

    rm results/$out/$id'_r.pdb'
done<$list

