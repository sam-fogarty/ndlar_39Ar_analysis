import numpy as np
import cvxpy as cp
from scipy.optimize import curve_fit

def multicomponent_fit(data, *MC):
    n = len(MC)
    A = np.zeros((len(data), n))
    
    for i, MC_sample in enumerate(MC):
        A[:, i] = MC_sample
    
    # construct the problem to solve
    x = cp.Variable(n)
    objective = cp.Minimize(cp.sum_squares(A @ x - data))
    constraints = [0 <= x, x >= 0]
    prob = cp.Problem(objective, constraints)
    
    # the optimal objective value is returned by prob.solve()
    result = prob.solve()
    
    # factors to multiply each component by
    values = x.value
    
    return values

def power_law(x, A, B):
    return A * x**(-B)

def fit_power_law(bincenters, data_bins, x_cut):
    bincenters_mask = bincenters >= x_cut
    bincenters_cut = bincenters[bincenters_mask]
    data_bins_cut = data_bins[bincenters_mask]
    
    popt, pcov = curve_fit(power_law, bincenters_cut, data_bins_cut)
    return power_law(bincenters, *popt), popt

import numpy as np
import cvxpy as cp
from scipy.optimize import curve_fit

def multicomponent_fit_with_powerlaw(data, *MC, bin_centers, A_init=1.0, B_init=1.0, tol=1e-6, max_iter=100):
    n = len(MC) + 1  # +1 for the power law component
    A = np.zeros((len(data), n))
    
    # Initialize A and B
    A_val = A_init
    B_val = B_init
    
    for iter_count in range(max_iter):
        # Update the components with the current A and B
        for i, MC_sample in enumerate(MC):
            A[:, i] = MC_sample
        A[:, -1] = A_val * bin_centers**(-B_val)  # power law component

        # Solve for the proportions using CVXPY
        x = cp.Variable(n)
        objective = cp.Minimize(cp.sum_squares(A @ x - data))
        constraints = [x >= 0]
        prob = cp.Problem(objective, constraints)
        result = prob.solve(solver=cp.SCS)
        x_val = x.value

        # Solve for A and B using curve_fit
        popt, pcov = curve_fit(power_law, bin_centers, data - np.dot(A[:, :-1], x_val[:-1]))
        A_val_new, B_val_new = popt

        # Check for convergence
        if np.abs(A_val_new - A_val) < tol and np.abs(B_val_new - B_val) < tol:
            break
        A_val, B_val = A_val_new, B_val_new

    # Return the proportions and the power law parameters
    return x_val, (A_val, B_val)
