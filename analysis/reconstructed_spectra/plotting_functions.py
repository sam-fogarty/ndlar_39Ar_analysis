import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import h5py
import scipy.stats

# constants
v_cm_sim = 288.28125
v_ref_sim = 1300.78125
v_pedestal_sim = 580 

# set ADC parameters for data and MC
v_cm_data = 288.28125
v_ref_data = 1300.78125

gain_data = 1/221
gain_sim = 0.004

def norm_hist(array, bins, range_start, range_end, norm, scale=1):
    ### make normalized histogram
    
    # get histogram data
    y,binEdges = np.histogram(np.array(array),bins=bins,range=(range_start,range_end))
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    menStd     = np.sqrt(y)
    y_norm = np.zeros_like(y,dtype='float64')
    y_norm_std = np.zeros_like(y,dtype='float64')
    
    # normalize bincontents
    bincontents_total = np.sum(y)
    for i in range(len(y)):
        y_uncert = ufloat(y[i],menStd[i])
        if norm == 'area':
            y_uncert = y_uncert/bincontents_total
        elif norm == 'max':
            y_uncert = y_uncert/np.max(y)
        else:
            y_uncert = y_uncert
        y_uncert *= scale
        y_norm[i] = y_uncert.nominal_value
        y_norm_std[i] = y_uncert.std_dev
    return bincenters, y_norm, y_norm_std

def get_hist_data(the_data, bins, data_type, calibrate=False, norm='none',binwidth=None,recomb_filename=None):
    ### set bin size and histogram range, correct for recombination using NEST, return histogram parameters
    # INPUT: `the_data` is a 1D numpy array of charge cluster charge values
    #        `bins` is the number of bins to use (this effectively sets the range, the binsize is constant w.r.t. bins)
    #        `data_type` is either `data` for real data or `MC` for simulation
    #        `calibrate` is either True or False, if True will use NEST to correct for recombination
    #        `norm` sets the normalization of histogram. Options are `area` and `max`, otherwise `None`
    #        `binwidth` is optional to specify binwidth. Otherwise will be 2*LSB by default
    #        `recomb_filename` is path to h5 file containing electron recombination values as a function of energy
    data = np.copy(the_data)
    
    # set parameters for data or MC for determining bin size
    vcm_mv = v_cm_data
    vref_mv = v_ref_data
    gain = gain_data
    
    if data_type == 'MC':
        vcm_mv = v_cm_sim
        vref_mv = v_ref_sim
        gain = gain_sim
    LSB = (vref_mv - vcm_mv)/256
    width = LSB / gain * 1e-3
    
    # add small +/- offset in MC for fix binning issues
    offset = np.zeros_like(data)
    if data_type == 'MC':
        MC_size = len(data)
        offset = scipy.stats.uniform.rvs(loc=0, scale=1, size=MC_size)*width*np.random.choice([-0.5,0.5], size=MC_size, p=[.5, .5])
        data += offset
    if binwidth is not None:
        width = binwidth
    
    range_start = width*0 # ke-
    range_end = width*bins
    
    # if converting to energy, use NEST model
    eV_per_e = 1
    if calibrate:
        eV_per_e = 23.6
        if recomb_filename == None:
            raise Exception("Calibrate is set to True, so must provide non-Null filename of h5 file with energies and recombination values.")
        else:
            recomb_file = h5py.File(recomb_filename)
            energies = np.array(recomb_file['NEST']['E_start'])
            recombination = np.array(recomb_file['NEST']['R'])
            charge_ke = energies / 23.6 * recombination
            data = data/np.interp(data, charge_ke, recombination)
    
    # get histogram parameters
    scale=1
    bincenters, y_norm, y_norm_std = norm_hist(data*eV_per_e, int(bins/2), range_start*eV_per_e,range_end*eV_per_e,norm,scale)
    return bincenters, y_norm, y_norm_std, data*eV_per_e

def plot_hist(bincenters, y_norm, y_norm_std, plots, color, linewidth, label):
    ### add spectrum histogram to matplotlib axes
    plots.step(bincenters, y_norm, linewidth=linewidth, color=color,where='mid',alpha=0.7, label=label)
    plots.errorbar(bincenters, y_norm, yerr=y_norm_std,color='k',fmt='o',markersize = 1)

    