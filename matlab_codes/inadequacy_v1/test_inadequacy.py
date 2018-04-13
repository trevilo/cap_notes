
import numpy as np
import matplotlib.pyplot as plt

# set parameters for both models
alpha = 0.33;
#alpha = 0.25;
alphaPrime = 1.0;
lamL  = np.pi*np.pi;
tauS  = 1.0/(4.0*np.pi*np.pi);
#tauS  = 1.0/(2.0*np.pi*np.pi);
lam2 = 4.0*np.pi*np.pi
sigma = 3.0;
beta  = 0.0002;
#beta  = 0.00005;
gam = 1.0;

# First, the step current case
Nsample = 64;

# Loading HF data
A = np.loadtxt("step_current/hf_lf_models/error_in_qoi.txt", delimiter=",");
tau       = A[:,0];
I         = A[:,1];
eps_exact = A[:,2]; 

# change first point
I[0] = 0.0
#I[0] = I[1]
eps_exact[0] = 0.0

dt = tau[1] - tau[0]
print "dt = ", dt
print "sqrt(dt) = ", np.sqrt(dt)
epsM = np.zeros((tau.shape[0], Nsample))
lamS = np.zeros((tau.shape[0], Nsample))

dW = np.reshape(np.sqrt(dt)*np.random.randn(tau.shape[0]*Nsample), [tau.shape[0], Nsample])

Ip = np.zeros(I.shape)

#beta  = 1.0e-4*np.random.randn(Nsample)
#beta  = 1.0e-4*np.random.rand(Nsample)
beta  = 1.5e-1*np.random.rand(Nsample)
gamma  = 2.0e-3*np.random.rand(Nsample)

lamFast = np.exp(1.0 + np.random.randn(Nsample))


for ii in range(1,tau.shape[0]):
    for jj in range(0,Nsample):
        #lamS[ii,jj] = lamS[ii-1,jj] - lamS[ii-1,jj]*(dt/tauS) + beta*np.abs(I[ii] - I[ii-1]) + sigma*lamS[ii-1,jj]*dW[ii-1,jj]
        #epsM[ii,jj] = epsM[ii-1,jj] - lamL*epsM[ii-1,jj]*dt - lamS[ii-1,jj]*epsM[ii-1,jj]*dt - alpha*(I[ii] - I[ii-1])
        #epsM[ii,jj] = epsM[ii-1,jj] - lamL*epsM[ii-1,jj]*dt - lamS[ii-1,jj]*epsM[ii-1,jj]*dt - alpha*(I[ii] - I[ii-1]) - alphaPrime*(I[ii] - I[ii-1])*np.sqrt(dt)

        #lamS[ii,jj] = lamS[ii-1,jj] - lam2*lamS[ii-1,jj]*dt - beta*(I[ii] - I[ii-1]) - 1000.0*epsM[ii-1,jj]**3*dt 
        #lamS[ii,jj] = lamS[ii-1,jj] - lam2*lamS[ii-1,jj]*dt - beta*(I[ii] - I[ii-1]) - 1000.0*epsM[ii-1,jj]**3*dt - 2.0*lamS[ii-1,jj] * dW[ii-1,jj]
        #lamS[ii,jj] = lamS[ii-1,jj] - lam2*lamS[ii-1,jj]*dt - beta*(I[ii] - I[ii-1])  - 0.25*(5.0*lamS[ii-1,jj] + 1000.0*epsM[ii-1,jj]**3) * dW[ii-1,jj] 
        #epsM[ii,jj] = epsM[ii-1,jj] - lamL*epsM[ii-1,jj]*dt - lamS[ii-1,jj]*dt - alpha*(I[ii] - I[ii-1])

        Ip[ii] = Ipp = rat = 0.0
        if (ii == 1 or ii == 2):
            Ip[ii]  = (I[ii+1] - I[ii-1])/(2.0*dt)
            Ipp = 0.0
            rat = np.abs(Ipp) #/np.maximum(np.abs(Ip[ii]),1e-8) 
        elif (ii < tau.shape[0]-1):
            Ip[ii]  = (I[ii+1] - I[ii-1])/(2.0*dt)
            Ipp = (I[ii+1] - 2.0*I[ii] + I[ii-1])/(dt*dt)
            rat = np.abs(Ipp) #/np.maximum(np.abs(Ip[ii]),1e-8) 
        else:
            Ip[ii]  = (I[ii] - I[ii-1])/(dt)
            Ipp = (I[ii] - 2.0*I[ii-1] + I[ii-2])/(dt*dt)
            rat = np.abs(Ipp) #/np.maximum(np.abs(Ip[ii]),1e-8) 

        print Ip[ii], Ipp, rat
        #lamS[ii,jj] = lamS[ii-1,jj] - 4.0*(lamL + lamS[ii-1,jj])*lamS[ii-1,jj]*dt + 0.01*rat*rat*dt
        #epsM[ii,jj] = epsM[ii-1,jj] - (lamL + lamS[ii-1,jj])*epsM[ii-1,jj]*dt - alpha*(I[ii] - I[ii-1])


        rat = np.minimum(rat, 0.5/dt)
        #print tau[ii], Ip[ii], Ipp, rat
        lamS[ii,jj] = lamS[ii-1,jj] - 4.0*(lamL + lamS[ii-1,jj])*lamS[ii-1,jj]*dt + 0.05*rat*rat*dt + dW[ii-1,jj]


        #print Ip[ii], Ipp, rat
        #lamS[ii,jj] = lamS[ii-1,jj] - (lamL + lamFast[jj])*lamS[ii-1,jj]*dt - beta[jj]*Ipp + 0.1*dW[ii-1,jj]
        ##lamS[ii,jj] = lamS[ii-1,jj] - (lamL)*lamS[ii-1,jj]*dt - beta[jj]*Ipp + 0.1*dW[ii-1,jj]
        #epsM[ii,jj] = epsM[ii-1,jj] - lamL*epsM[ii-1,jj]*dt - alpha*(I[ii] - I[ii-1]) + lamS[ii-1,jj]*dt
        epsM[ii,jj] = epsM[ii-1,jj] - (lamL + lamS[ii-1,jj])*epsM[ii-1,jj]*dt - alpha*(I[ii] - I[ii-1]) + beta[jj]*Ip[ii]*dt - gamma[jj]*Ipp*dt





mean_epsM = np.mean(epsM, axis=1)
std_epsM = np.std(epsM, axis=1)
min_epsM = np.amin(epsM, axis=1)
max_epsM = np.amax(epsM, axis=1)

plt.figure()
#plt.plot(tau, epsM[:,0], 'b-', linewidth=2, label="Model")
plt.plot(tau, mean_epsM, 'b-', linewidth=2, label="Model mean")
plt.plot(tau, min_epsM, 'b-.', linewidth=2, label="Model bounds")
plt.plot(tau, max_epsM, 'b-.', linewidth=2)
for jj in range(0,Nsample):
    plt.plot(tau, epsM[:,jj], 'k-', linewidth=1)

plt.plot(tau, eps_exact, 'r--', linewidth=2, label="Exact")
plt.xlabel(r"$\tau$", fontsize=16)
plt.ylabel(r"$\epsilon$", fontsize=16)
plt.ylim(-0.3, 0.3)
plt.legend()
plt.grid()
plt.savefig("epsilon_step_current.pdf", bbox_inches="tight")
plt.close()

plt.figure()
plt.plot(tau, lamS[:,0], 'b-', linewidth=2, label="Model")
for jj in range(0,Nsample):
    plt.plot(tau, lamS[:,jj], 'k-', linewidth=1)

plt.xlabel(r"$\tau$", fontsize=16)
plt.ylabel(r"$\lambda_S$", fontsize=16)
plt.legend()
plt.grid()
plt.savefig("lamS_step_current.pdf", bbox_inches="tight")
plt.close()

rhsH = np.zeros(tau.shape)
lamE = np.zeros(tau.shape)
for ii in range(1,tau.shape[0]-1):
    dEps = (eps_exact[ii+1] - eps_exact[ii-1])/(2.0*dt)
    dI = (I[ii+1] - I[ii-1])/(2.0*dt)
    rhsH[ii] = dEps + lamL*eps_exact[ii] + alpha*(I[ii] - I[ii-1])
    lamE[ii] = rhsH[ii]/eps_exact[ii]
    
plt.figure()
plt.plot(tau[:-1], rhsH[:-1], 'b-')
plt.plot(tau[:-1], lamS[:-1], 'r-')
#plt.plot(tau[:-1], Ip[:-1], 'g--')
plt.savefig('lam_exact_step.pdf', bbox_inches='tight')


# Next, the sinusoidal current case
Nsample = 64;
gamma = np.logspace(-4, 0, Nsample)

# Loading HF data
A = np.loadtxt("sinusoid_current/hf_lf_models/error_in_qoi_sin.txt", delimiter=",");
tau       = A[:,0];
I         = A[:,1];
eps_exact = A[:,2]; 

dt = tau[1] - tau[0]
epsM = np.zeros((tau.shape[0], Nsample))
lamS = np.zeros((tau.shape[0], Nsample))

dW = np.reshape(np.sqrt(dt)*np.random.randn(tau.shape[0]*Nsample), [tau.shape[0], Nsample])

##beta  = np.logspace(-4, -1, Nsample)
#beta  = np.linspace(-1e-4, 1e-4, Nsample)
#beta  = 1.0e-4*np.random.randn(Nsample)
#beta  = np.linspace(0,1,Nsample)
#gamma  = np.linspace(0,1,Nsample)

#beta  = 2.0e-1*np.random.rand(Nsample)
#gamma  = 2.0e-3*np.random.rand(Nsample)

#beta  = 1.5e-1*np.random.rand(Nsample)
#beta  = 0.1*np.ones(Nsample)
#gamma  = 4.0e-3*np.random.rand(Nsample)
beta  = 1.0e-1*np.random.randn(Nsample)
gamma  = 2.0e-1*np.random.rand(Nsample)


lamFast = np.exp(1.0 + np.random.randn(Nsample))
##beta  = np.logspace(-5, -3, Nsample)

#beta = 0.00005
#lamS[0,:] = 0.0
lamS[0,:] = -0.5

for ii in range(1,tau.shape[0]):
    for jj in range(0,Nsample):
        #lamS[ii,jj] = lamS[ii-1,jj] - lamS[ii-1,jj]*(dt/tauS) + beta*np.abs(I[ii] - I[ii-1]) + sigma*lamS[ii-1,jj]*dW[ii-1,jj]
        #epsM[ii,jj] = epsM[ii-1,jj] - lamL*epsM[ii-1,jj]*dt - lamS[ii-1,jj]*epsM[ii-1,jj]*dt - alpha*(I[ii] - I[ii-1]) - alphaPrime*(I[ii] - I[ii-1])*np.sqrt(dt)

        #lamS[ii,jj] = lamS[ii-1,jj] - lam2*lamS[ii-1,jj]*dt - beta*(I[ii] - I[ii-1]) - 1000.0*epsM[ii-1,jj]**3*dt 
        #lamS[ii,jj] = lamS[ii-1,jj] - lam2*lamS[ii-1,jj]*dt - beta*(I[ii] - I[ii-1])  - 2.0*lamS[ii-1,jj] * dW[ii-1,jj] 
        #lamS[ii,jj] = lamS[ii-1,jj] - lam2*lamS[ii-1,jj]*dt - beta*(I[ii] - I[ii-1])  - 0.25*(5.0*lamS[ii-1,jj] + 1000.0*epsM[ii-1,jj]**3) * dW[ii-1,jj] 
        #epsM[ii,jj] = epsM[ii-1,jj] - lamL*epsM[ii-1,jj]*dt - lamS[ii-1,jj]*dt - alpha*(I[ii] - I[ii-1])

        #Ip[ii] = Ipp = rat = 0.0
        #if (ii < tau.shape[0]-1):
        #Ip[ii]  = (I[ii+1] - I[ii-1])/(2.0*dt)
        #Ipp = (I[ii+1] - 2.0*I[ii] + I[ii-1])/(dt*dt)
        #rat = np.abs(Ipp)/np.maximum(np.abs(Ip[ii]),1e-8) 
        ##rat = np.abs(Ipp)/(4.0)
        #else:
        #Ip[ii]  = (I[ii] - I[ii-1])/(dt)
        #Ipp = (I[ii] - 2.0*I[ii-1] + I[ii-2])/(dt*dt)
        #rat = np.abs(Ipp)/np.maximum(np.abs(Ip[ii]),1e-8) 
        ##rat = np.abs(Ipp)/(4.0) 
        #
        #rat = np.minimum(rat, 0.5/dt)
        #print tau[ii], Ip[ii], Ipp, rat
        ##lamS[ii,jj] = lamS[ii-1,jj] - 4.0*(lamL + lamS[ii-1,jj])*lamS[ii-1,jj]*dt + 0.00001*rat*rat*dt
        #lamS[ii,jj] = lamS[ii-1,jj] - 4.0*(lamL + lamS[ii-1,jj])*lamS[ii-1,jj]*dt + gamma[jj]*rat*rat*dt
        ##lamS[ii,jj] = lamS[ii-1,jj] - (lamL + lamS[ii-1,jj])*lamS[ii-1,jj]*dt + 0.1*rat*rat*dt
        ##lamS[ii,jj] = lamS[ii-1,jj] - lamS[ii-1,jj]*lamS[ii-1,jj]*dt + 0.5*rat*rat*dt
        ##lamS[ii,jj] = lamS[ii-1,jj] - lamS[ii-1,jj]*dt + gam*rat*rat*dt
        ##lamS[ii,jj] = lamS[ii-1,jj] - lamS[ii-1,jj]*lamS[ii-1,jj]*dt
        #print "ii = ", ii, ", lamS = ", lamS[ii,jj]
        ##epsM[ii,jj] = epsM[ii-1,jj] - (lamL + lamS[ii-1,jj])*epsM[ii-1,jj]*dt - alpha*(I[ii] - I[ii-1])
        ##epsM[ii,jj] = epsM[ii-1,jj] - (lamL + lamS[ii-1,jj])*epsM[ii-1,jj]*dt - alpha*(I[ii] - I[ii-1])
        ##epsM[ii,jj] = epsM[ii-1,jj] - (lamL + lamS[ii-1,jj])*epsM[ii-1,jj]*dt - 0.2*(I[ii] - I[ii-1])
        #epsM[ii,jj] = epsM[ii-1,jj] - (lamL + lamS[ii-1,jj])*epsM[ii-1,jj]*dt - alpha*(I[ii] - I[ii-1])
        ##epsM[ii,jj] = epsM[ii-1,jj] - (lamL+20.0)*epsM[ii-1,jj]*dt - alpha*(I[ii] - I[ii-1])
        
        Ip[ii] = Ipp = rat = 0.0
        if (ii < tau.shape[0]-1):
            Ip[ii]  = (I[ii+1] - I[ii-1])/(2.0*dt)
            Ipp = (I[ii+1] - 2.0*I[ii] + I[ii-1])/(dt*dt)
            rat = np.abs(Ipp)#/np.maximum(np.abs(Ip[ii]),1e-8) 
        else:
            Ip[ii]  = (I[ii] - I[ii-1])/(dt)
            Ipp = (I[ii] - 2.0*I[ii-1] + I[ii-2])/(dt*dt)
            rat = np.abs(Ipp)#/np.maximum(np.abs(Ip[ii]),1e-8) 

        print Ip[ii], Ipp, rat
        #lamS[ii,jj] = lamS[ii-1,jj] - 4.0*(lamL + lamS[ii-1,jj])*lamS[ii-1,jj]*dt + 0.01*rat*rat*dt
        #lamS[ii,jj] = lamS[ii-1,jj] - 4.0*(lamL + lamS[ii-1,jj])*lamS[ii-1,jj]*dt + 0.1*rat*rat*dt
        lamS[ii,jj] = lamS[ii-1,jj] - 4.0*(lamS[ii-1,jj])*lamS[ii-1,jj]*dt + 0.1*rat*rat*dt
        #epsM[ii,jj] = epsM[ii-1,jj] - (lamL + beta[jj]*lamS[ii-1,jj])*epsM[ii-1,jj]*dt - alpha*(I[ii] - I[ii-1])


        rat = np.minimum(rat, 0.5/dt)
        #print tau[ii], Ip[ii], Ipp, rat
        #lamS[ii,jj] = lamS[ii-1,jj] - 4.0*(lamL + lamS[ii-1,jj])*lamS[ii-1,jj]*dt + 0.01*rat*rat*dt
        lamS[ii,jj] = lamS[ii-1,jj] - 4.0*(lamL + lamS[ii-1,jj])*lamS[ii-1,jj]*dt + 0.01*rat*rat*dt + dW[ii-1,jj]

        epsM[ii,jj] = epsM[ii-1,jj] - (lamL + lamS[ii-1,jj])*epsM[ii-1,jj]*dt - alpha*(I[ii] - I[ii-1]) + beta[jj]*Ip[ii]*dt - gamma[jj]*Ipp*dt


mean_epsM = np.mean(epsM, axis=1)
std_epsM = np.std(epsM, axis=1)
min_epsM = np.amin(epsM, axis=1)
max_epsM = np.amax(epsM, axis=1)

plt.figure()
#plt.plot(tau, epsM[:,0], 'b-', linewidth=2, label="Model")
#plt.plot(tau, epsM[:,0], 'k-', linewidth=1)
for jj in range(0,Nsample):
    plt.plot(tau, epsM[:,jj], 'k-', linewidth=1)
#plt.plot(tau, epsM[:,20], 'k-', linewidth=1)
#plt.plot(tau, epsM[:,30], 'k-', linewidth=1)
#plt.plot(tau, epsM[:,40], 'k-', linewidth=1)
plt.plot(tau, mean_epsM, 'b-', linewidth=2, label="Model mean")
plt.plot(tau, min_epsM, 'b-.', linewidth=2, label="Model bounds")
plt.plot(tau, max_epsM, 'b-.', linewidth=2)
plt.plot(tau, eps_exact, 'r--', linewidth=2, label="Exact")
#plt.plot(tau, I, 'g--', linewidth=2, label="Current")
#plt.plot(tau, Ip, 'g--', linewidth=2, label="Current")
plt.xlabel(r"$\tau$", fontsize=16)
plt.ylabel(r"$\epsilon$", fontsize=16)
plt.legend(loc=4)
plt.grid()
plt.savefig("epsilon_sine_current.pdf", bbox_inches="tight")
plt.close()

plt.figure()
plt.plot(tau, lamS[:,0], 'b-', linewidth=2, label="Model")
for jj in range(0,Nsample):
    plt.plot(tau, lamS[:,jj], 'k-', linewidth=1)
#plt.plot(tau, Ip, 'g--', linewidth=2, label="Current")
plt.xlabel(r"$\tau$", fontsize=16)
plt.ylabel(r"$\lambda_S$", fontsize=16)
plt.legend(loc=0)
plt.grid()
plt.savefig("lamS_sine_current.pdf", bbox_inches="tight")
plt.close()


rhsH = np.zeros(tau.shape)
lamE = np.zeros(tau.shape)
for ii in range(1,tau.shape[0]-1):
    dEps = (eps_exact[ii+1] - eps_exact[ii-1])/(2.0*dt)
    dI = (I[ii+1] - I[ii-1])/(2.0*dt)
    Ipp = (I[ii+1] - -2.0*I[ii] + I[ii-1])/(dt*dt)
    rhsH[ii] = dEps + lamL*eps_exact[ii] + 0.7*alpha*dI
    lamE[ii] = rhsH[ii]/Ipp
    
plt.figure()
plt.plot(tau[:-1], rhsH[:-1], 'b-')
#plt.plot(tau[:-1], lamE[:-1], 'b-')
plt.plot(tau[:-1], lamS[:-1], 'g--')
plt.savefig('lam_exact_sine.pdf', bbox_inches='tight')
