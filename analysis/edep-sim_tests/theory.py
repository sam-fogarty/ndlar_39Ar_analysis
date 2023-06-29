## script to calculate the track lengths for particles according to theory
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from theory_fxns import *

# integrand for finding full mean range of particle
def FullRange(E):
    return 1/(dEdx_mean(E))

def dEdxANDRR(E_i, N):
    #E_i = 0.565 # MeV
    E_end = 0.01 # small, ideally zero but might cause problems
    tracklength = abs(integrate.quad(FullRange, E_i, E_end)[0]) # cm

    #N = 500
    dx = tracklength/N
    energies = np.linspace(E_end, E_i, N)
    residualRange = np.arange(0.0, tracklength, dx)
    dEdx = np.array([dEdx_mean(E) for E in energies])
    return residualRange, dEdx, energies
