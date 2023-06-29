from theory_consts import *
import numpy as np

def gamma(E_KE):
    return 1 + E_KE/M

def beta(E_KE):
    return np.sqrt(1 - 1/(gamma(E_KE)**2))

def L2(E_KE):
    term_sum = 0
    y = z/(137*beta(E_KE))
    for l in range(1,steps):
        term_sum += 1/(l**3 /y*y + l)
    return -1*z*z*term_sum

def Wmax(E_KE):
    #return 2*m_e*beta(E_KE)**2 * gamma(E_KE)**2 /(2 + 2*gamma(E_KE))
    return m_e*(gamma(E_KE)-1)

def dEdx_mean_mine(E_KE):
    # my calculation based on PDG
    # for delta effect correction
    x = np.log10(beta(E_KE)*gamma(E_KE))
    dcorr = 0
    if x >= x1:
        dcorr = 2*np.log(10)*x + C
    elif x >= x0 and x < x1:
        dcorr = 2*np.log(10)*x + C + a*(x1 - x)**k
    elif x < x0:
        dcorr = 0
    
    # for L1 correction. Make sure to define F function.
    y = (137*beta(E_KE))**2 / Z
    #L1 = np.sqrt(2) * F(b/(y**(1/2)) / (Z**(1/2) * y**(3/2))
    L1 = np.sqrt(2) * F / (Z**(1/2) * y**(3/2))
    L1=0
    # calculate bracketed terms of bethe-bloche eqn, eqn 33.5 in pdg
    bracket = np.log(2*m_e*beta(E_KE)**2 * gamma(E_KE)**2 * Wmax(E_KE) /(I*I))/2 - beta(E_KE)**2 - dcorr/2
    bracket = bracket + L2(E_KE) + L1
    
    dEdx_mean = -1*bracket*K*z*z*Z/A * 1/beta(E_KE)**2
    return abs(dEdx_mean * rho)

def dEdx_mean(E_KE):
    # see eqn 34.23 in PDG
    b = beta(E_KE)
    g = gamma(E_KE)
    # for density effect correction
    x = np.log10(b*g)
    dcorr = 0
    if x >= x1:
        dcorr = 2*np.log(10)*x + C
    elif x >= x0 and x < x1:
        dcorr = 2*np.log(10)*x + C + a*(x1 - x)**k
    elif x < x0:
        dcorr = 0
        
    logterm = np.log(m_e*b**2 * g**2 * (m_e*(g-1)/2)/I**2)
    bracket = logterm + (1 - b**2) - ((2.0*g - 1.0)/g**2)*np.log(2) + (1/8)*((g-1)/g)**2 - dcorr
    return bracket*rho*0.5*z*z*K*(Z/A)*b**(-2)