This repository contains all scripts needed to run the bioinformatics analysis for
the paper:

Required input files should be downloaded from the specified locations:

run_analysis.sh details how all scripts and commands are used throughout the analysis.
It is set up so that it could run everything automatically provided all the necessary
dependencies and input files are in place.

In particular, to complete the miRNA prediction and annotation part of the pipeline,
miRCat, miRDeep and MapMi output files are required. Copies of these outputs
are provided in the directory ./miRNAs. These can be used as input to the script
process_annotations.sh to merge and process all miRNA annotations.
