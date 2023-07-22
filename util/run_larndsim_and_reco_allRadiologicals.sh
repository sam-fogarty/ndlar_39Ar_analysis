#!/bin/bash

# define constants
input_folder=/sdf/group/neutrino/sfogarty/ND_prototype_files/MC/module-0/radiologicals_cachedLArNEST/edep-sim
nevents=10000
filenum=3
output_folder=/sdf/group/neutrino/sfogarty/ND_prototype_files/MC/module-0/radiologicals_cachedLArNEST/larnd-sim
output_folder_new=/sdf/group/neutrino/sfogarty/ND_prototype_files/MC/module-0/radiologicals_cachedLArNEST/reco
original_dir=$(pwd)
larndsim_dir=/sdf/home/s/sfogarty/Desktop/RadDecay/larnd-sim_NEST/larnd-sim_cachedLArNEST_radiologicals/cli
log_file=$original_dir/loop_runtime_$filenum.log

# define the list of isotope and decay combinations
declare -A isotopes
isotopes["39Ar"]="betas"
isotopes["85Kr"]="betas gammas"
isotopes["60Co"]="betas gammas"
isotopes["40K"]="betas gammas"
isotopes["232Th"]="betas gammas alphas"
isotopes["238U"]="betas gammas alphas"

# record the start time
start_time=$SECONDS
echo "Script started at: $(date)" >> $log_file

cd $larndsim_dir

# loop over isotopes
for isotope in "${!isotopes[@]}"
do
    IFS=" "
    decays=(${isotopes[$isotope]})

    # loop over decay modes
    for decay in "${decays[@]}"
    do
	# record the start time of the loop
        loop_start_time=$SECONDS
        
	# generate filenames
        input_filename="edep_${isotope}_${decay}_${nevents}_${filenum}.h5"
        output_filename="larndsim_${isotope}_${decay}_${nevents}_${filenum}.h5"
        output_filename_new="larndsim_${isotope}_${decay}_${nevents}_${filenum}_clusters.h5"

        # run simulate_pixels.py
        python3 simulate_pixels.py \
            --input_filename=$input_folder/$input_filename \
            --output_filename=$output_folder/$output_filename \
            --detector_properties=../larndsim/detector_properties/module0.yaml \
            --simulation_properties=../larndsim/simulation_properties/singles_sim.yaml \
            --pixel_layout=../larndsim/pixel_layouts/multi_tile_layout-2.4.16.yaml \
            --response_file=../larndsim/bin/response_44.npy \
            --pixel_thresholds_file=module0-measured_thresholds-6ke.npz
        # change to reco directory and run reco.py
        cd /sdf/home/s/sfogarty/Desktop/RadDecay/ndlar_hitfinder/reco
        python3 reco.py module0_MC.py $output_folder/$output_filename $output_folder_new/$output_filename_new
	
	cd $original_dir
	# record the end time of the loop and calculate the loop runtime
        loop_end_time=$SECONDS
        echo "Loop for $isotope $decay finished at: $(date). Runtime: $((loop_end_time-loop_start_time)) seconds." >> $log_file

        # return to original directory
        cd $larndsim_dir
    done

done

cd $original_dir

# record the end time and calculate the total runtime
end_time=$SECONDS
echo "Script ended at: $(date). Total runtime: $((end_time-start_time)) seconds." >> $log_file
