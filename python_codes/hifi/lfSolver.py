import numpy as np
from scipy.integrate import cumtrapz


def solve(tau, xi, I, gamma):

    eta_bar = cumtrapz(I, x=tau, initial=tau[0]*I[0])
    
    etaLF = np.zeros((xi.shape[0], tau.shape[0]))
    for jj in range(0,xi.shape[0]):
        for ii in range(0,tau.shape[0]):
            etaLF[jj,ii] = ((I[ii]/2)*xi[jj]**2 - 
                            ((I[ii]*gamma)/(1.0+gamma))*xi[jj] + 
                            eta_bar[ii] - I[ii]/6.0 + 
                            (I[ii]*gamma)/(2*I[ii]+2*gamma))

    return etaLF
