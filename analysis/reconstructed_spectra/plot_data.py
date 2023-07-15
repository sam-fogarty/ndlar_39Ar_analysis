#!/usr/bin/env python

import fire
import numpy as np
import h5py
import matplotlib.pyplot as plt
import plotting_functions

def main(filepath):
    # open data h5 file and get clusters
    f = h5py.File(filepath, 'r')
    charge_data = f['clusters']['q']
    
    recomb_filename = '/sdf/home/s/sfogarty/Desktop/RadDecay/LArNDLE/sim/larnd-sim/NEST/NEST_R-values_efield500_1keV_to_10000keV_1keV-stepsize_1000-events.h5'
    fig, axes = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, figsize=(8,5))
    
    # make data histogram
    do_calibration = False
    normalization = 'None'
    bincenters_data, y_norm_data, y_norm_std_data, charge_calib_data = \
        plotting_functions.get_hist_data(charge_data, bins=75, data_type='data', \
        calibrate=do_calibration, norm=normalization,recomb_filename=recomb_filename)
    plotting_functions.plot_hist(bincenters_data, y_norm_data, y_norm_std_data, axes, color='b', linewidth=1, label='data')
    
    if do_calibration:
        axes.set_xlabel('Energy [keV]')
    else:
        axes.set_xlabel('Charge [ke-]')
    
    axes.set_title('Reconstructed 39Ar beta decay spectrum (Module-0)')
    axes.legend()
    plt.show()
    
if __name__ == "__main__":
    fire.Fire(main)