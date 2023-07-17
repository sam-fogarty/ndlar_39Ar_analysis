import numpy as np
import fire
import larnestpy
import h5py
import matplotlib.pyplot as plt
from tqdm import tqdm

def main(start_keV, end_keV, step_size):
    # Calculate recombination factors (1 - recombination probability) from LArNEST for electron recoils
    # for a range of energies and store in a file.
    # Inputs:
    #     `start_keV`: first energy in keV to calculate recombination factor
    #     `end_keV`: last energy in keV to calculate recombination factor
    #     `step_size`: step size in keV
    # The code creates a text file and an h5py file containing the data, and a pdf of the energies vs R.
    energy_vals = np.arange(start_keV,end_keV+step_size, step_size, dtype='float32')
    calc = larnestpy.LArNEST()
    
    R = np.zeros_like(energy_vals)
    
    # calculate recombination factor for each energy
    efield = 500 # V/cm
    dx = 1
    density = 1.393 # g/cm^3
    for i, E in enumerate(energy_vals):
        result = calc.full_calculation(larnestpy.LArInteraction.ER, E, dx, efield, density, False)
        R[i] = 1 - result.yields.RecombinationProbability
    
    # save recombination values to a text file for larnd-sim
    filepath = f'NEST_R-values_efield_{efield}_{start_keV}_keV_to_{end_keV}_keV_{step_size}_keV_stepsize'
    print('Writing text file ', filepath+'.txt')
    f = open(filepath+'.txt',"w+")
    for i, R in enumerate(R):
        E = energy_vals[i]
        f.write(str(E) + ',' + str(R))
        f.write('\n')
    f.close()
    
    # load in text file
    data = np.loadtxt(filepath+'.txt', dtype='float', delimiter=',')
    energies = data[:,0]
    recombination = data[:,1]
    
    # make plot and save pdf
    plt.plot(energies, recombination)
    plt.xlabel('Energy [keV]')
    plt.ylabel('Recombination Factor R')
    plt.savefig(filepath+'.pdf')
    
    filepath_h5 = filepath+'.h5'
    print('Writing h5 file ', filepath_h5)
    # save recombination values and energies in h5 file
    dtype = np.dtype([('E_start', 'f4'),('R', 'f4')])
    datapoints = np.empty(len(energies), dtype=dtype)
    for i in range(len(energies)):
        datapoints[i]['E_start'] = energies[i]
        datapoints[i]['R'] = recombination[i]

    with h5py.File(filepath_h5, 'w') as f:
        f.create_dataset("NEST", data=datapoints)
    
if __name__ == "__main__":
    fire.Fire(main)