%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Script to generate HF data
%%%          for constant current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all


%% Parameters
kappa = 0.0195174;   %--- sec/m
sigma = 52.1;   %--- sec/m
gamma = kappa/sigma;
L = 50e-6;  %--- m
C = 0.03134;  %--- F/m2
a = 4.19956e7/C;   %--- m
V0 = 1.25;   %--- volt
Iunscaled = 200;   %--- Amp/m^2
Ls = 25e-6;  %--- m
kappa_s = 0.0311627;   %--- sec/m

Nt = 5000;  %--- number of time steps
Nx = 700;  %--- number of space steps 

%% Domain
t = linspace(0.001,5 , Nt);  %--- sec
x = linspace(0,L , Nx);

% convert to nondimesional values
tau = (kappa*sigma/(kappa+sigma))*(1/(a*C*L^2)) .* t;
xi = x ./ L;
Nxi = Nx ; Ntau = Nt;

%% Applied Current
t_conv = a*C*L^2*(kappa+sigma)/(kappa*sigma);

Iamp = (Iunscaled*(L/V0)*((kappa+sigma)/(sigma*kappa)));
Ins = Iamp*square(6*pi*tau');
I = Ins';
% I = smooth(Ins,'moving');

%%========================================================
%% Computing Vcell
%%========================================================


%% Solution of LF model

% compute eta avg
eta_bar(1) = tau(1)*I(1);
for i = 2 : size(tau,2)
    eta_bar(i) = trapz(tau(1:i),I(1:i));
end

for i = 1 : length(tau) % time
    for j = 1 : length(xi) % space
        etaLF(j,i) = (I(i)/2)*xi(j)^2 - ((I(i)*gamma)/(1+gamma))*xi(j) + eta_bar(i) - I(i)/6 + (I(i)*gamma)/(2*I(i)+2*gamma);
    end
end


etaLF_xi0 = etaLF(1,:);
etaLF_xi1 = etaLF(end,:);

V_elecLF = ((1+2*gamma)/(1+gamma))*etaLF_xi1 - ...
    (gamma/(1+gamma))*etaLF_xi0 - (gamma/(1+gamma)^2)*I;


%% Exact Solution of HF model 
C_HF = -(gamma/(1+gamma)^2);
V_elecHF = C_HF*I;
for i = 1:length(I)
    V_elecHF(i) = V_elecHF(i) + tIntegral(tau(i)); % int K_hf*I
end

B = (Ls/kappa_s)*(kappa*sigma/(L*(kappa+sigma)));


%% QoI
V_cellHF = 1 - 0.5*B*I - V_elecHF;
V_cellLF = 1 - 0.5*B*I - V_elecLF;

figure
plot(tau,V_cellHF,'b','LineWidth',3); hold on
plot(tau,V_cellLF,'r','LineWidth',3);
legend('V_{HF}','V_{LF}')
xlabel('\tau')
ylabel('$$QoI = V^{\rm cell}$$','Interpreter','LaTex')
prop_plots


%% Error in QoI

eps_exact =  V_cellHF - V_cellLF;

figure
plot(tau,eps_exact,'LineWidth',3)
xlabel('\tau')
ylabel('\epsilon')
legend('HF data')
prop_plots

%%========================================================
%% Save in data file
M = [tau' I' eps_exact'];
csvwrite('error_in_qoi_stepI.txt',M);
