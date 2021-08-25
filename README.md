# BindingSite_docking

```
Steps to reproduce data:
1 - Clone this directory


2 - Enter the data folder and unpack all zipped files.

3 - Download Benchmark4 dataset (set4) from this URL:
    http://dockground.bioinformatics.ku.edu/unbound/request.php
    Move it in the data directory and fully expand it.

4 - Download the Gramm output from: [figshare]
    fully expand it and move all the content of the resulting folder in BindingSite_docking/results/

From step 5 you'll need several python modules, look at the Singularity section of this README for more details.

5 - From the directory BindingSite_docking run:

    python3 src/process_benchmark4.py -c list_dimers -i data/benchmark4/pdb/ -o data/processed_b4
    
    OPTIONAL - if you want to recreate the formatted_labels:
    5A - From the directory BindingSite_docking run:
         
         python3 src/format_dataset.py -s [path to BindingSite_docking folder on your system]
         
         You may modify this script to format your own Binding Sites predictions and test them into the pipeline.
       
       
6 - To select docking poses with binding sites predictions run from the directory BindingSite_docking:
    
    python3 src/protocol_gramm.py \
        -g results/gramm_output/[id]_u1-[id]_u2.res \
        -i1 data/formatted_labels/[id]_u1.site \
        -i2 data/formatted_labels/[id]_u2.site \
        -s1 data/processed_b4/[id]_u1.pdb \
        -s2 data/processed_b4/[id]_u2.pdb \
        -c [col] -n [#CPU] -o [output_dir]
    
    substituting [id] with a pdb 4 letter code available in processed_b4, [col] with the column of the formatted_labels
    file you want to use for the scoring (counted starting from 0), [#CPU] with the number of CPU you want to use to 
    parallelize and [output_dir] with the path and prefix where you want the output files to be saved.
    You can take a look at the src/run_gramm.sh script to have a sample on how to launch this script.
    
7 - Take a look at https://github.com/bjornwallner/DockQ to obtain the DockQ program which is used in this study to 
    assess docking quality. You can find the pre-computed results in the results folder, named summary_dockq_[col], 
    where [col] is a number referring to the column in formatted_labels files used in protocol_gramm.py to select 
    the related dockings. Notice that these are also column indexes and they go from 1 to 18 (first column in the files
    indexed 0 is amino acid solvent accessibility).
    
8 - Open the jupyter notebook (src/plot_notebook.ipynb) and set correctly the paths in the first cell in order to run the 
    data analysis.
    
SINGULARITY ENVIRONMENT
In order to easily provide all the necessary modules (mainly Tensorflow and Biopython) two singularity recipe files are 
available in the src folder (singularity_def_GPU and singularity_def_CPU in case no GPU is available). A list of all 
dependencies can be found in these files. To use these singularity definition files to create an environment and run 
commands check: https://sylabs.io/guides/3.0/user-guide/

        


```
