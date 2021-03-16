# BindingSite_docking

```
Steps to reproduce data:
1 - Clone this directory

2 - Enter the folder and unpack the data directory

3 - Download Benchmark4 dataset (set4) in the data directory
    from this URL: http://dockground.bioinformatics.ku.edu/unbound/request.php
    and unpack it
    
4 - run:

    python3 process_benchmark4.py -c list_dimers -i data/benchmark4/pdb/ -o data/processed_b4

5 . run:

    while read id; do for n in {1..2}; do python3 format_dataset.py -c ${id}_u${n}.pdb -s data/processed_b4/ -r data/naccess/ -o data/formatted_labels/ -isp data/ISPRED_predictions/ -dyn data/dynJET2_predictions/ -spp data/SPPIDER_predictions/; done; done<list_dimers

NOTES
CONSURFDB RUN (hopefully): https://consurfdb.tau.ac.il/scripts/waitBatch.php?runNumber=1615331348

removed 1rpq_u2 and 4cu4_u2 cuz too short for SPPIDER
removed 2o3b_u2 for failure in dynJET2

dynJET2 ran for 3 iterations in mode 4

```
