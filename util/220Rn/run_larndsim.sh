#!/bin/bash

# define constants
input_folder=/sdf/group/neutrino/sfogarty/ND_prototype_files/MC/module-0/edep-sim_tests
nevents=10000
filenum=1
output_folder=/sdf/group/neutrino/sfogarty/ND_prototype_files/MC/module-0/larnd-sim_tests
original_dir=$(pwd)
larndsim_dir=/sdf/home/s/sfogarty/Desktop/RadDecay/larnd-sim/cli

isotope='39Ar'
decay='betas'
output_file_descriptor='randSeed12345'
# record the start time
start_time=$SECONDS
echo "Script started at: $(date)"

cd $larndsim_dir

# record the start time of the loop
loop_start_time=$SECONDS

# generate filenames
#input_filename="edep_${isotope}_${decay}_${nevents}_events_1.h5"
#output_filename="larndsim_${isotope}_${decay}_${nevents}_events_${output_file_descriptor}_${filenum}.h5"
input_filename="MiniRun4_1E19_RHC.convert2h5.00000.EDEPSIM.h5"
output_filename="MiniRun4_1E19_RHC.convert2h5.00000.EDEPSIM_larndsim_1.h5"
# run simulate_pixels.py
python3 simulate_pixels.py \
    --input_filename=$input_folder/$input_filename \
    --output_filename=$output_folder/$output_filename \
    --detector_properties=../larndsim/detector_properties/2x2.yaml \
    --simulation_properties=../larndsim/simulation_properties/2x2_NuMI_sim.yaml \
    --pixel_layout=../larndsim/pixel_layouts/multi_tile_layout-2.4.16.yaml \
    --response_file=../larndsim/bin/response_44.npy #\
    #--rand_seed=12345
#--pixel_thresholds_file=module0-measured_thresholds-6ke.npz
      
#--detector_properties=../larndsim/detector_properties/module0.yaml \
#--simulation_properties=../larndsim/simulation_properties/singles_sim.yaml \
cd $original_dir

# record the end time and calculate the total runtime
end_time=$SECONDS
echo "Script ended at: $(date). Total runtime: $((end_time-start_time)) seconds." 
