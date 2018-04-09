# Script to plot inadequecy representation
import numpy as np
import scipy.signal as spsig
import matplotlib.pyplot as plt

# files don't exist on git repo
#eps_exact = csvread('error_cycl.txt');
#error_exact_Velec = csvread('error_cycl.txt');
#V_elecHF = csvread('V_elecHF.txt');
#V_elecLF = csvread('V_elecLF.txt');

# Inadequacy Parameters
Nsamp = 17;
lambda_mean = 15;
#c = 15;
c = 30;
#alpha = 0.32;
alpha = 0.32;
#zeta = 0.0002;
#zeta = 0.0002*0.5;
zeta = 0.00016*0.25

#beta = 15;
beta = 0.2;


# Parameters
kappa = 0.0195174;   #--- sec/m
sigma = 52.1;   #--- sec/m
gamma = kappa/sigma;
L = 50e-6;  #--- m
C = 0.03134;  #--- F/m2
a = 4.19956e7/C;   #--- m
V0 = 1.25;   #--- volt
Iunscaled = 200;   #--- Amp/m^2
Ls = 25e-6;  #--- m
kappa_s = 0.0311627;   #--- sec/m
Time = 4; #--- sec

#Nt = 800;  #--- number of time steps
#Nt = 1600;  #--- number of time steps
Nt = 3200;  #--- number of time steps

# Initial conditions
Epsilon = np.zeros((Nt, Nsamp))
Epsilon[0,:] = 0.02*np.random.randn(Nsamp)+0.14;

Lambda = np.zeros((Nt, Nsamp))
Lambda[0, :] = 20.0;



# Domain
t = np.linspace(0,Time , num=Nt);  #--- sec

# convert to nondimesional values
tau = (kappa*sigma/(kappa+sigma))*(1/(a*C*L*L)) * t;
Ntau = Nt;

# Applied Current
t_conv = a*C*L*L*(kappa+sigma)/(kappa*sigma);

Iamp = (Iunscaled*(L/V0)*((kappa+sigma)/(sigma*kappa)));
Ins = Iamp*spsig.square(6*np.pi*tau);
I = Ins;

dtau = tau[2] - tau[1];

# Inadequecy evolution
dW = np.sqrt(dtau)*np.reshape(np.random.randn(Nsamp*Ntau), (Ntau, Nsamp));   #increments

print Lambda.shape
print Epsilon.shape

dIdt = np.zeros(Lambda.shape[0])
force = np.zeros(Lambda.shape[0])

for i in range(0,Ntau-2):
    for j in range(0,Nsamp):
        #Lambda[i+1,j]  = Lambda[i,j]  - c*(Lambda[i,j]-lambda_mean)*dtau + beta*dW[i,j];
        #Epsilon[i+1,j] = Epsilon[i,j] -    Lambda[i,j]*Epsilon[i,j]*dtau + alpha*(I[i+1]-I[i]);
        #Epsilon[i+1,j] = Epsilon[i,j] -    Lambda[i,j]*Epsilon[i,j]*dtau + 0.5*alpha*(I[i+2]-I[i]) + zeta*(I[i+2] - 2.0*I[i+1] + I[i])/dtau;

        dIdt[i+1] = 0.5*(I[i+2] - I[i])/dtau
        force[i+1] = beta*np.sqrt(np.abs(dIdt[i+1]))*np.exp(100*dW[i,j])
        
        Lambda[i+1,j]  = Lambda[i,j]  - c*(Lambda[i,j]-lambda_mean)*dtau + force[i+1]
        #Lambda[i+1,j]  = Lambda[i,j]  - c*(Lambda[i,j]-lambda_mean)*dtau
        Epsilon[i+1,j] = Epsilon[i,j] -    Lambda[i,j]*Epsilon[i,j]*dtau + 0.5*alpha*(I[i+2]-I[i]);



# figure
# plot(tau,lambda,'--r','LineWidth',3)
# xlabel('\tau'); ylabel('\lambda');
# prop_plots
# 
# figure
# plot(tau,error_exact_Velec,'b','LineWidth',3); hold on
# plot(tau,-Epsilon,'--','LineWidth',2)
# xlabel('\tau'); ylabel('Epsilon');
# prop_plots
#%
#figure
#plot(tau,V_elecHF,'r','LineWidth',3); hold on
#plot(tau,V_elecLF,'--b','LineWidth',3);
#xlabel('Normalized Time \tau'); ylabel('QoI: V^{cell}');
#legend ('V_{HF}','V_{LF}')
#prop_plots
#daspect([1 2 1])

# confidence plot
Meps = np.mean(-Epsilon, axis=1);
Leps = np.amin(-Epsilon, axis=1);
Ueps = np.amax(-Epsilon, axis=1);

plt.figure()
plt.fill_between(tau, Leps, Ueps, color="grey", alpha=0.5)
plt.plot(tau, Meps, 'b-', linewidth=2, label="mean")
plt.plot(tau, -Epsilon[:,0], 'g-', linewidth=1, label="Sample0")
plt.plot(tau, -Epsilon[:,1], 'g-', linewidth=1, label="Sample1")
plt.savefig('epsilon_results.pdf', bbox_inches='tight')



# confidence plot
Mlam = np.mean(Lambda, axis=1);
Llam = np.amin(Lambda, axis=1);
Ulam = np.amax(Lambda, axis=1);

plt.figure()
plt.fill_between(tau, Llam, Ulam, color="grey", alpha=0.5)
plt.plot(tau, Mlam, 'b-', linewidth=2, label="mean")
plt.plot(tau, Lambda[:,0], 'g-', linewidth=1, label="Sample0")
plt.plot(tau, Lambda[:,1], 'g-', linewidth=1, label="Sample1")
plt.savefig('lambda_results.pdf', bbox_inches='tight')



# confidence plot
plt.figure()
plt.plot(tau, dIdt, 'b-', linewidth=2, label="mean")
plt.savefig('dIdt_results.pdf', bbox_inches='tight')

plt.figure()
plt.plot(tau, force, 'b-', linewidth=2, label="mean")
plt.savefig('force.pdf', bbox_inches='tight')


#%% plot QoI correction
#for i = 1 : Nsamp
#    Vlf_corrected(i,:) = V_elecLF(:)-Epsilon(:,i);
#end
#Mv = mean(Vlf_corrected);
#Lv = min(Vlf_corrected);
#Uv = max(Vlf_corrected);
#
#figure
#p = plot(tau,Lv,tau, Uv,tau,Mv);
#YLIM = get(gca,'YLim'); delete(p);
#a1 = area(tau,Uv,min(YLIM)); 
#hold on
#set(a1,'LineStyle','none');     set(a1,'FaceColor',[0.9 0.9 0.9]);
#a2 = area(tau,Lv,min(YLIM)); 
#set(a2,'LineStyle','none');     set(a2,'FaceColor',[1 1 1]);
#p_mean = plot(tau,Mv,'--b','LineWidth',2);
#p_exact = plot(tau,V_elecHF,'r','LineWidth',3);
#xlabel('Normalized Time \tau'); ylabel('QoI: V^{cell}');
#legend([p_exact p_mean a1],{'V_{HF}','V_{LF}+error','95% confidence boundaries'})
#prop_plots
#daspect([1 2 1])
