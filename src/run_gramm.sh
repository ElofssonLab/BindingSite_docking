list=$1
out=$2
col=$3

folder=`pwd`

if [[ ! -f $out ]]
then
    mkdir results/$out
fi

while read id; do
    echo 'Processing ' $id '...'
    singularity exec --bind $folder:$folder --nv $folder/BSenv.sif python $folder/src/protocol_gramm.py -s1 $folder/data/processed_b4/$id'_u1.pdb' -s2 $folder/data/processed_b4/$id'_u2.pdb' -o $folder/results/$out/$id -g $folder/data/gramm-res/$id'_u1-'$id'_u2.res' -i1 $folder/data/formatted_labels/$id'_u1.site' -i2 $folder/data/formatted_labels/$id'_u2.site' -c $col;
    for n in {1..10}; do 
	cat $folder/data/processed_b4/$id'_u1.pdb' | grep -v 'END' > $folder/results/$out/$id'_'$n'-tmp'; 
	echo '' >> $folder/results/$out/$id'_'$n'-tmp';  
	cat $folder/results/$out/$id'_'$n'.pdb' >> $folder/results/$out/$id'_'$n'-tmp'; 
	mv $folder/results/$out/$id'_'$n'-tmp' $folder/results/$out/$id'_'$n'.pdb'; 
    done;
done<$list
