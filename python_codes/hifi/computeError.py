import numpy as np
import hfSolver as hf
import lfSolver as lf

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.signal import square

# Set up constant parameters
    
# Parameters
kappa = 0.0195174;   #--- sec/m
sigma = 52.1;        #--- sec/m
gamma = kappa/sigma;
L = 50e-6;           #--- m
C = 0.03134;         #--- F/m2
a = 4.19956e7/C;     #--- m
V0 = 1.25;           #--- volt
Iunscaled = 200;     #--- Amp/m^2
Ls = 25e-6;          #--- m
kappa_s = 0.0311627; #--- sec/m
Time = 4;            #--- sec
    
Nt = 1024;           #--- number of time steps
Nx = 256;            #--- number of space steps


## Domain
t = np.linspace(0,Time , Nt);  #--- sec
x = np.linspace(0,L , Nx);

# convert to nondimesional values
tau = (kappa*sigma/(kappa+sigma))*(1.0/(a*C*L*L)) * t;
xi = x / L;
Nxi = Nx
Ntau = Nt

## Applied Current
Iamp = (Iunscaled*(L/V0)*((kappa+sigma)/(sigma*kappa)));


## Different cases

#========================================================
# Square wave case
print "Running step wave case..."
I = Iamp*square(6.0*np.pi*tau);

# Prep the BCs
alpha = -I*(gamma/(1+gamma));
beta  =  I*(1/(1+gamma));

# IC
eta0 = np.zeros(Nx)

# hf solve
xg, tg, etaHF = hf.solve(tau[0], tau[-1], xi[0], xi[-1], Ntau, Nxi, alpha, beta, eta0, plot=False)
qoiHF = ((1+2*gamma)/(1+gamma))*etaHF[-1,:] - (gamma/(1+gamma))*etaHF[0,:] - (gamma/(1+gamma)**2)*I;

# lf solve
etaLF = lf.solve(tau, xi, I, gamma)
qoiLF = ((1+2*gamma)/(1+gamma))*etaLF[-1,:] - (gamma/(1+gamma))*etaLF[0,:] - (gamma/(1+gamma)**2)*I;

# error
err = qoiHF - qoiLF

W = np.zeros((Nt,3))
W[:,0] = tau
W[:,1] = I
W[:,2] = err
np.savetxt("step_qoi_error.txt", W)

plt.plot(tau, err, 'b-', linewidth=2)
plt.savefig('qoi_error_step.pdf', bbox_inches='tight')



#========================================================
# Sine 1
print "Running sine wave case (1)..."
I = Iamp*np.sin(5.0*np.pi*tau);

# Prep the BCs
alpha = -I*(gamma/(1+gamma));
beta  =  I*(1/(1+gamma));

# IC
eta0 = np.zeros(Nx)

# hf solve
xg, tg, etaHF = hf.solve(tau[0], tau[-1], xi[0], xi[-1], Ntau, Nxi, alpha, beta, eta0, plot=False)
qoiHF = ((1+2*gamma)/(1+gamma))*etaHF[-1,:] - (gamma/(1+gamma))*etaHF[0,:] - (gamma/(1+gamma)**2)*I;

# lf solve
etaLF = lf.solve(tau, xi, I, gamma)
qoiLF = ((1+2*gamma)/(1+gamma))*etaLF[-1,:] - (gamma/(1+gamma))*etaLF[0,:] - (gamma/(1+gamma)**2)*I;

# error
err = qoiHF - qoiLF

W = np.zeros((Nt,3))
W[:,0] = tau
W[:,1] = I
W[:,2] = err
np.savetxt("sin1_qoi_error.txt", W)

plt.figure()
plt.plot(tau, err, 'b-', linewidth=2)
plt.savefig('qoi_error_sin1.pdf', bbox_inches='tight')


#========================================================
# Sine 1
print "Running sine wave case (2)..."
I = Iamp*np.sin(10.0*np.pi*tau);

# Prep the BCs
alpha = -I*(gamma/(1+gamma));
beta  =  I*(1/(1+gamma));

# IC
eta0 = np.zeros(Nx)

# hf solve
xg, tg, etaHF = hf.solve(tau[0], tau[-1], xi[0], xi[-1], Ntau, Nxi, alpha, beta, eta0, plot=False)
qoiHF = ((1+2*gamma)/(1+gamma))*etaHF[-1,:] - (gamma/(1+gamma))*etaHF[0,:] - (gamma/(1+gamma)**2)*I;

# lf solve
etaLF = lf.solve(tau, xi, I, gamma)
qoiLF = ((1+2*gamma)/(1+gamma))*etaLF[-1,:] - (gamma/(1+gamma))*etaLF[0,:] - (gamma/(1+gamma)**2)*I;

# error
err = qoiHF - qoiLF

W = np.zeros((Nt,3))
W[:,0] = tau
W[:,1] = I
W[:,2] = err
np.savetxt("sin2_qoi_error.txt", W)

plt.figure()
plt.plot(tau, err, 'b-', linewidth=2)
plt.savefig('qoi_error_sin2.pdf', bbox_inches='tight')

