#!/bin/bash

# define constants
nevents=10000
physics_list=QGSP_BERT_LIV
output_folder=/sdf/group/neutrino/sfogarty/ND_prototype_files/MC/module-0/radiologicals_afterNESTfix/edep-sim/
filenum=5
dumpTree_dir=~/Desktop/RadDecay/larnd-sim_NEST/larnd-sim_cache/cli

# define the list of isotope, isotope_full, decay combinations
declare -A isotopes
isotopes["39Ar"]="Argon_39,betas"
isotopes["85Kr"]="Krypton_85,betas gammas"
isotopes["60Co"]="Cobalt_60,betas gammas"
isotopes["40K"]="Potassium_40,betas gammas"
isotopes["232Th"]="Thorium_232,betas gammas alphas"
isotopes["238U"]="Uranium_238,betas gammas alphas"

# loop over isotopes
for isotope in "${!isotopes[@]}"
do
	    IFS=","
	        read -a full_and_decays <<< "${isotopes[$isotope]}"
    isotope_full="${full_and_decays[0]}"

    IFS=" "
    decays=(${full_and_decays[1]})

    # loop over decay modes
    for decay in "${decays[@]}"
    do
        # generate filenames
        edep_filename="edep_${isotope}_${decay}_${nevents}_${filenum}.root"
        macro_file="${isotope_full}_${decay}.mac"

        # execute edep-sim
        edep-sim -e $nevents -p $physics_list -o $output_folder/${edep_filename} ../macros/${macro_file}

        # convert the output file
        python3 $dumpTree_dir/dumpTree.py $output_folder/${edep_filename} $output_folder/edep_${isotope}_${decay}_${nevents}_${filenum}.h5
    done
done

