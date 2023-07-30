#!/bin/bash
### bash script for running ndlar_flow LRS reconstruction

original_dir=$(pwd)

# change to your own directories
ndlar_flow_dir=/global/u1/s/sfogarty/ndlar_flow
output_folder=/global/cfs/cdirs/dune/users/sfogarty

# change which files are used and which module is used
data_folder=/global/cfs/cdirs/dune/www/data/Module1/LRS/SingleModule_Jan22
ADC_name_1=0cd913fb
ADC_name_2=0cd93db0
file_timestamp=20220208_014759

file_ending=data
LRS_file_1=${ADC_name_1}_${file_timestamp}.${file_ending}
LRS_file_2=${ADC_name_2}_${file_timestamp}.${file_ending}
light_workflow_yamls=yamls/module0_flow/workflows/light
light_reco_yamls=yamls/module0_flow/reco/light
output_file=${output_folder}/${ADC_name_1}_${ADC_name_2}_${file_timestamp}.h5

# you can tweak these files as needed before copying over
cp light_yamls/light_event_reconstruction.yaml ${ndlar_flow_dir}/${light_workflow_yamls}
cp light_yamls/WaveformHitFinder.yaml ${ndlar_flow_dir}/${light_reco_yamls}

cd ${ndlar_flow_dir}
pip install .

# run light event building adc64 step for both files
h5flow -c ${light_workflow_yamls}/light_event_building_adc64.yaml \
    -i ${data_folder}/${LRS_file_1} \
    -o ${output_file}
    
h5flow -c ${light_workflow_yamls}/light_event_building_adc64.yaml \
    -i ${data_folder}/${LRS_file_2} \
    -o ${output_file}

# run light event reconstruction step
h5flow -c ${light_workflow_yamls}/light_event_reconstruction.yaml \
    -i ${output_file} \
    -o ${output_file}
