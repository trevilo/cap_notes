import numpy as np
import matplotlib.pyplot as plt

# load step data
data = np.loadtxt("./step_current/inadequecy/error_in_qoi.txt", delimiter=",")
tau = data[:,0]
I   = data[:,1]
err = data[:,2]

Nt = err.shape[0]
dtau = tau[1] - tau[0]

Dt = np.zeros((Nt,Nt))

for ii in range(1,Nt-1):
    Dt[ii,ii+1] = 0.5/dtau;
    Dt[ii,ii-1] = -0.5/dtau;

Dt[0,1] = 1.0/dtau;
Dt[0,0] = -1.0/dtau;

Dt[-1,-1] = 1.0/dtau;
Dt[-1,-2] = -1.0/dtau;

alfa = 0.333;
derrdt = np.inner(Dt,err) - alfa*np.inner(Dt,I)

lam = np.pi*np.pi

#beta = (-derrdt/err - lam)/(err*err+1e-16)
beta = (-derrdt/err - lam)

#plt.plot(tau, beta, 'b-', label='step')
              

# load step data
data = np.loadtxt("./sinusoid_current/hf_lf_models/error_in_qoi_sin.txt", delimiter=",")
tau = data[:,0]
I   = data[:,1]
err = data[:,2]

Nt = err.shape[0]
dtau = tau[1] - tau[0]

Dt = np.zeros((Nt,Nt))

for ii in range(1,Nt-1):
    Dt[ii,ii+1] = 0.5/dtau;
    Dt[ii,ii-1] = -0.5/dtau;

Dt[0,1] = 1.0/dtau;
Dt[0,0] = -1.0/dtau;

Dt[-1,-1] = 1.0/dtau;
Dt[-1,-2] = -1.0/dtau;

alfa = 0.333;
#alfa = 0.25;
#alfa = 0.2;
dIdt = np.inner(Dt,I)
derrdt = np.inner(Dt,err) + alfa*np.inner(Dt,I);
#derrdt = np.inner(Dt,err)

lam = np.pi*np.pi

beta1 = (derrdt)
beta2 = (derrdt + lam*err)
beta3 = (derrdt + lam*err + 200*err*err*err)

print "norm beta2 = ", np.linalg.norm(beta2)
print "norm beta3 = ", np.linalg.norm(beta3)

plt.plot(tau, dIdt, 'k-', label='step')
plt.plot(tau, beta1, 'r-', label='step')
plt.plot(tau, beta2, 'b-', label='step')
plt.plot(tau, beta3, 'g-', label='step')

plt.grid()
#plt.ylim(-3000,3000)              
plt.savefig('test_beta.pdf', bbox_inches='tight')

