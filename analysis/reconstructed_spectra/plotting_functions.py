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
#v_cm_data = 284.2773
#v_ref_data = 1282.7148

gain_data = 1/221
#gain_data = 0.004
#gain_sim = 0.004
gain_sim = 1/221
def make_hist(array, bins, range_start, range_end):
    ### make histogram of charge
    
    # get histogram data
    bin_contents,binEdges = np.histogram(np.array(array),bins=bins,range=(range_start,range_end))
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    error      = np.sqrt(bin_contents)
    
    return bincenters, bin_contents, error

def get_hist_data(the_data, bins, data_type, calibrate=False, binwidth=None, recomb_filename=None):
    ### set bin size and histogram range, correct for recombination using NEST, return histogram parameters
    # INPUT: `the_data` is a 1D numpy array of charge cluster charge values (mV)
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
    
    data = data/gain * 1e-3
    
    if data_type == 'MC':
        vcm_mv = v_cm_sim
        vref_mv = v_ref_sim
        gain = gain_sim
    LSB = (vref_mv - vcm_mv)/256
    if calibrate:
        width = LSB / gain * 1e-3
    else:
        width = LSB
    
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
    bin_centers, bin_contents, bin_error = make_hist(data*eV_per_e, int(bins/2), range_start*eV_per_e,range_end*eV_per_e)
    return bin_centers, bin_contents, bin_error

def plot_hist(bin_centers, bin_contents, bin_error, plots, color, linewidth, label, linestyle=None, norm=None,dont_plot_errorbars=False):
    ### add spectrum histogram to matplotlib axes
    
    if norm == 'area':
        total_bin_contents = np.sum(bin_contents)
        bin_contents = bin_contents / total_bin_contents
        bin_error = bin_error / total_bin_contents

    plots.step(bin_centers, bin_contents, linewidth=linewidth, color=color,linestyle=linestyle, where='mid',alpha=0.7, label=label)
    if not dont_plot_errorbars:
        plots.errorbar(bin_centers, bin_contents, yerr=bin_error,color='k',fmt='o',markersize = 1)
    
def get_charge_MC(nFiles_dict, folders_MC, filename_ending_MC, nbins, do_calibration, recomb_filename,disable_alphas=False, disable_gammas=False, disable_betas=False):
    # Isotope ratios
    isotopes_ratios_betas_gammas = { 
        '85Kr': 224.76, # beta/gamma ratio
        '60Co': 0.5,
        '40K': 8.46
    }
    
    isotopes_ratios_betas_alphas_gammas_alphas = {
        '232Th': [0.649, 0.45], # betas/alphas , gammas/alphas ratios
        '238U': [0.751, 0.999]
    }
    
    # Initialize dictionaries
    charge_dict = {}
    hist_data_dict = {}
    
    # Loop over isotopes
    for iso_decay, nFiles in nFiles_dict.items():
        iso, decay = iso_decay.split('_')
        folder = folders_MC[iso_decay]
        ending = filename_ending_MC[iso_decay]
        # Loop over files
        for i in range(1, nFiles+1):
            f = h5py.File(folder + f'larndsim_{iso}_{decay}_10000_{i}_{ending}.h5', 'r')
            charge_temp = f['clusters']['q']
            if i == 1:
                charge_dict[iso_decay] = charge_temp
            else:
                charge_dict[iso_decay] = np.concatenate((charge_dict[iso_decay], charge_temp))
        
        # Call function to get histogram data
        bin_centers, bin_contents, bin_error = \
            get_hist_data(charge_dict[iso_decay], bins=nbins, data_type='MC', \
            calibrate=do_calibration, recomb_filename=recomb_filename)
        
        hist_data_dict[iso_decay] = {
            'bin_centers': bin_centers,
            'bin_contents': bin_contents,
            'bin_error': bin_error
        }
    
    # Combine y_norm and y_norm_std for isotopes that have betas and gammas
    for iso in isotopes_ratios_betas_gammas.keys():
        R = isotopes_ratios_betas_gammas[iso]
        x_1 = R * (np.sum(hist_data_dict[iso+'_gammas']['bin_contents']) / np.sum(hist_data_dict[iso+'_betas']['bin_contents']))
        x_2 = 1
        
        if disable_gammas:
            x_2 = 0
        if disable_betas:
            x_1 = 0
        hist_data_dict[iso] = {
            'bin_centers': hist_data_dict[iso+'_betas']['bin_centers'],
            'bin_contents': hist_data_dict[iso+'_betas']['bin_contents']*x_1 + hist_data_dict[iso+'_gammas']['bin_contents']*x_2,
            'bin_error': np.sqrt((hist_data_dict[iso+'_betas']['bin_error']*x_1)**2 + (hist_data_dict[iso+'_gammas']['bin_error']*x_2)**2)
        }
    for iso in isotopes_ratios_betas_alphas_gammas_alphas.keys():
        R = isotopes_ratios_betas_alphas_gammas_alphas[iso]
        x_1 = R[0] * (np.sum(hist_data_dict[iso+'_alphas']['bin_contents']) / np.sum(hist_data_dict[iso+'_betas']['bin_contents']))
        x_2 = R[1] * (np.sum(hist_data_dict[iso+'_alphas']['bin_contents']) / np.sum(hist_data_dict[iso+'_gammas']['bin_contents']))
        x_3 = 1
        if disable_betas:
            x_1 = 0
        if disable_gammas:
            x_2 = 0
        if disable_alphas:
            x_3 = 0
        hist_data_dict[iso] = {
            'bin_centers': hist_data_dict[iso+'_betas']['bin_centers'],
            'bin_contents': hist_data_dict[iso+'_betas']['bin_contents']*x_1 + hist_data_dict[iso+'_gammas']['bin_contents']*x_2 + \
                (hist_data_dict[iso+'_alphas']['bin_contents']*x_3),
            'bin_error': np.sqrt((hist_data_dict[iso+'_betas']['bin_error']*x_1)**2 + (hist_data_dict[iso+'_gammas']['bin_error']*x_2)**2 \
                                + (hist_data_dict[iso+'_alphas']['bin_error']*x_3)**2)
        }
    return charge_dict, hist_data_dict

def plot_isotopes(hist_data_dict, axes, colors, norm=None, linewidth=2, do_not_plot_list=None):    
    # Loop over isotopes
    for iso_decay, color in colors.items():
        # Get histogram data
        bin_centers = hist_data_dict[iso_decay]['bin_centers']
        bin_contents = hist_data_dict[iso_decay]['bin_contents']
        bin_error = hist_data_dict[iso_decay]['bin_error']
        
        if len(iso_decay.split('_')) > 1:
            label = iso_decay.split('_')[0]
        else:
            label = iso_decay
        if label not in do_not_plot_list:
            # Call function to plot histogram
            plot_hist(bin_centers, bin_contents, bin_error, axes, color, linewidth, label, norm=norm)
    
