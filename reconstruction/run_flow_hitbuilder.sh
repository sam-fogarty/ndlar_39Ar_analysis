#!/bin/bash
### bash script for running ndlar_flow charge hit building

ndlar_flow_dir=/global/u1/s/sfogarty/ndlar_flow
original_dir=$(pwd)
#data_folder=/global/cfs/cdirs/dune/www/data/Module1/TPC12/dataRuns/packetData
data_folder=/global/cfs/cdirs/dune/www/data/Module0/TPC1+2/dataRuns/packetData
output_folder=/global/cfs/cdirs/dune/users/sfogarty
#input_packets_filename=packet_2022_02_08_01_47_59_CET
input_packets_filename=datalog_2021_04_04_20_59_11_CEST
file_ending=.h5
file_descriptor=hits
charge_workflows_dir=yamls/module0_flow/workflows/charge
charge_building_yaml=charge_event_building.yaml
hit_building_yaml=charge_hit_building.yaml

cp ./${hit_building_yaml} ${ndlar_flow_dir}/${charge_workflows_dir}
source /global/homes/s/sfogarty/flow.venv/bin/activate

# run ndlar_flow event building step
cd $ndlar_flow_dir
pip install .
h5flow -c ${charge_workflows_dir}/${charge_building_yaml} -i ${data_folder}/${input_packets_filename}${file_ending} -o ${output_folder}/${input_packets_filename}_${file_descriptor}${file_ending}

# run_ndlar_flow hit building
h5flow -c ${charge_workflows_dir}/${hit_building_yaml} -i ${output_folder}/${input_packets_filename}_${file_descriptor}${file_ending} -o ${output_folder}/${input_packets_filename}_${file_descriptor}${file_ending}

cd ${original_dir}