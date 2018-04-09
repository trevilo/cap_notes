# Compute and plot the kernel, just for examination
import numpy as np
import matplotlib.pyplot as plt

#function integrand = k(t,tau)
kappa = 0.0195174;   #--- sec/m
sigma = 52.1;   #--- sec/m
gamma = kappa/sigma;

tt = np.logspace(-4, 1, num=1e3)

# My kernel
K = (2.0*gamma*gamma + (1.0 + 2.0*gamma - gamma*gamma))*np.ones(tt.shape)/((1.0+gamma)**2)

Agam = 2.0*gamma*gamma/((1.0+gamma)**2)
Bgam = (1.0 + 2.0*gamma - gamma*gamma)/((1.0+gamma)**2)

for n in range(1,100):
    K = K + 2.0*(Agam*(-1)**n + Bgam)*np.exp(-np.pi*np.pi*n*n*tt)


integrand = np.ones(tt.shape)
for n in range(1,100):
    integrand = integrand + 2*(2*gamma/(1+gamma)**2*(-1)**n+(1+gamma**2)/(1+gamma)**2)*np.exp(-n**2*np.pi**2*(tt))



plt.loglog(tt, K, 'b-', linewidth=2, label='Todd')
plt.loglog(tt, integrand, 'r--', linewidth=2, label='Danial')
plt.grid()
plt.xlabel(r'$t$', fontsize=18)
plt.ylabel(r'$K$', fontsize=18)
plt.legend()
plt.savefig('kernel.pdf', bbox_inches='tight')

# Danial's original matlab code
#integrand = 1;
#for n=1:100
#    integrand = integrand + 2*(2*gamma/(1+gamma)^2*(-1)^n+(1+gamma^2)/(1+gamma)^2)...
#        *exp(-n^2*pi^2*(tau-t));
#end
#end

