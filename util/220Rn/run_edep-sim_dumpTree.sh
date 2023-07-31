#!/bin/bash

# define constants
nevents=20000
nevents_name=20k

physics_list=QGSP_BERT_LIV
output_folder=/sdf/group/neutrino/sfogarty/ND_prototype_files/MC/2x2
dumpTree_dir=~/Desktop/RadDecay/larnd-sim/cli
isotope='212Po'
isotope_full='Polonium_212'
decay='alphas'

# generate filenames
edep_filename="edep_${isotope}_${decay}_${nevents_name}.root"
macro_file="${isotope_full}_${decay}.mac"

# execute edep-sim
edep-sim -e $nevents -p $physics_list -o $output_folder/${edep_filename} ../macros/${macro_file}

# convert the output file
python3 $dumpTree_dir/dumpTree.py $output_folder/${edep_filename} $output_folder/edep_${isotope}_${decay}_${nevents_name}.h5
