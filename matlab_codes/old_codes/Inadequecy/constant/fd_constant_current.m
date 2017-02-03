%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script to provide data for calibrating 
%%%         model inadequecy
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

Nt = 700;  %--- number of time steps
Nx = 700;  %--- number of space steps 

%% Domain
t = linspace(0,3 , Nt);  %--- sec
x = linspace(0,L , Nx);

% convert to nondimesional values
tau = (kappa*sigma/(kappa+sigma))*(1/(a*C*L^2)) .* t;
xi = x ./ L;
Nxi = Nx ; Ntau = Nt;

%% Applied Current
t_conv = a*C*L^2*(kappa+sigma)/(kappa*sigma);

I = (Iunscaled*(L/V0)*((kappa+sigma)/(sigma*kappa))) * ones(1,Nt);

% figure(1)
% plot(tau,I);
% axis square
% xlabel('\tau'); ylabel('I');
% prop_plots


%%========================================================
%% Solution of LF model
eta_bar = tau.*I;

for i = 1 : size(tau,2) % time
    for j = 1 : size(xi,2) % space
        etaLF(j,i) = (I(i)/2)*xi(j)^2 - ((I(i)*gamma)/(1+gamma))*xi(j) + eta_bar(i) - I(i)/6 + (I(i)*gamma)/(2*I(i)+2*gamma);
    end
end
% 
% figure
% surf(tau, xi, etaLF,'EdgeColor','none')
% axis square
% xlabel('\tau'); ylabel('\xi'); zlabel('{\eta}_{LF}');
% prop_plots
% 
% figure
% plot(xi,etaLF(:,Nt*0.25),'--','LineWidth',3); hold on
% plot(xi,etaLF(:,Nt*0.5),'LineWidth',3);
% plot(xi,etaLF(:,Nt*0.75),'LineWidth',3);
% plot(xi,etaLF(:,Nt*1),'-o','LineWidth',3);
% legend('\tau=0.25\tau_{total}', '\tau=0.5\tau_{total}', '\tau=0.75\tau_{total}','\tau=\tau_{total}')
% axis square
% xlabel('\xi'); ylabel('{\eta}_{LF}');
% prop_plots

%%========================================================
%% Solution of HF model
% 
% % computing Neumann BCs:
%-- alpha at xi=0
alpha = -I.*(gamma/(1+gamma));
%-- beta at xi=1
beta  =  I.*(1/(1+gamma));

% solving HF pde+bc using finite difference
etaHF = HFm_evolution(tau(1),tau(end),xi(1),xi(end), Ntau,Nxi, zeros(Nxi,Ntau), alpha, beta);

% figure
% surf(tau, xi, etaHF,'EdgeColor','none')
% axis square
% xlabel('\tau'); ylabel('\xi'); zlabel('{\eta}_{LF}');
% prop_plots
% 
% figure
% plot(xi,etaHF(:,Nt*0.25),'--','LineWidth',3); hold on
% plot(xi,etaHF(:,Nt*0.5),'LineWidth',3);
% plot(xi,etaHF(:,Nt*0.75),'LineWidth',3);
% plot(xi,etaHF(:,Nt*1),'-o','LineWidth',3);
% legend('\tau=0.25\tau_{total}', '\tau=0.5\tau_{total}', '\tau=0.75\tau_{total}','\tau=\tau_{total}')
% axis square
% xlabel('\xi'); ylabel('{\eta}_{HF}');
% prop_plots
% 


%%========================================================
%% Computing Vcell
etaHF_xi0 = etaHF(1,:);
etaHF_xi1 = etaHF(end,:);
etaLF_xi0 = etaLF(1,:);
etaLF_xi1 = etaLF(end,:);

V_elecHF = ((1+2*gamma)/(1+gamma))*etaHF_xi1 - (gamma/(1+gamma))*etaHF_xi0 - (gamma/(1+gamma)^2)*I;
V_elecLF = ((1+2*gamma)/(1+gamma))*etaLF_xi1 - (gamma/(1+gamma))*etaLF_xi0 - (gamma/(1+gamma)^2)*I;

B = (Ls/kappa_s)*(kappa*sigma/(L*(kappa+sigma)));

%%========================================================
%% Error in Velectrod

error_exact_Velec = V_elecHF - V_elecLF; 

figure
plot(tau,abs(error_exact_Velec),'b','LineWidth',3); hold on
xlabel('\tau'); ylabel('Error V_{electrod}');
prop_plots
