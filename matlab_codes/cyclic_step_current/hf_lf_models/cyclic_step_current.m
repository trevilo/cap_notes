%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script to solve the simple example
%%%     for model inadequecy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc


%% Parameters
kappa = 0.0195174;   %--- sec/m
sigma = 52.1;        %--- sec/m
gamma = kappa/sigma;
L = 50e-6;           %--- m
C = 0.03134;         %--- F/m2
a = 4.19956e7/C;     %--- m
V0 = 1.25;           %--- volt
Iunscaled = 200;     %--- Amp/m^2
Ls = 25e-6;          %--- m
kappa_s = 0.0311627; %--- sec/m
Time = 4;            %--- sec

Nt = 2000;            %--- number of time steps
Nx = 2000;            %--- number of space steps


%% Domain
t = linspace(0,Time , Nt);  %--- sec
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
%I = smooth(Ins,'moving');

figure
plot(I,'LineWidth',3);
xlabel('\tau'); ylabel('I^*');
prop_plots


%%========================================================
%% Solution of LF model

% compute eta avg
eta_bar(1) = tau(1)*I(1);
for i = 2 : size(tau,2)
    eta_bar(i) = trapz(tau(1:i),I(1:i));
end
figure
plot(tau,eta_bar,'r','LineWidth',3);
xlabel('\tau'); ylabel('\eta_{avg}');
prop_plots


for i = 1 : size(tau,2) % time
    for j = 1 : size(xi,2) % space
        etaLF(j,i) = (I(i)/2)*xi(j)^2 - ((I(i)*gamma)/(1+gamma))*xi(j) + eta_bar(i) - I(i)/6 + (I(i)*gamma)/(2*I(i)+2*gamma);
    end
end

figure
surf(tau, xi, etaLF,'EdgeColor','none')
xlabel('\tau'); ylabel('\xi'); zlabel('{\eta}_{LF}');
prop_plots

figure
plot(xi,etaLF(:,Nt*0.25),'--','LineWidth',3); hold on
plot(xi,etaLF(:,Nt*0.5),'LineWidth',3); plot(xi,etaLF(:,Nt*0.75),'LineWidth',3);
plot(xi,etaLF(:,Nt*1),'-o','LineWidth',3);
legend('\tau=0.25\tau_{total}', '\tau=0.5\tau_{total}', '\tau=0.75\tau_{total}','\tau=\tau_{total}')
xlabel('\xi'); ylabel('{\eta}_{LF}');
prop_plots

figure
plot(tau,etaLF(Nx*0.25,:),'--','LineWidth',3); hold on
plot(tau,etaLF(Nx*0.5,:),'LineWidth',3); plot(tau,etaLF(Nx*0.75,:),'LineWidth',3);
plot(tau,etaLF(Nx*1,:),'-o','LineWidth',3);
legend('\xi=0.25\xi_{total}', '\xi=0.5\xi_{total}', '\xi=0.75\xi_{total}','\xi=\xi_{total}')
xlabel('\tau'); ylabel('{\eta}_{LF}');
prop_plots


%%========================================================
%% Solution of HF model
%
% Neumann BCs:
%-- alpha at xi=0
alpha = -I.*(gamma/(1+gamma));
%-- beta at xi=1
beta  =  I.*(1/(1+gamma));

% IC:
eta0(:,1) = 0;


% solving HF pde+bc using finite difference
etaHF = fdm(tau(1),tau(end),xi(1),xi(end), Ntau,Nxi, alpha, beta, eta0);

figure
surf(tau, xi, etaHF,'EdgeColor','none')
xlabel('\tau'); ylabel('\xi'); zlabel('{\eta}_{LF}');
prop_plots


figure
plot(xi,etaHF(:,Nt*0.25),'--','LineWidth',3); hold on
plot(xi,etaHF(:,Nt*0.5),'LineWidth',3); plot(xi,etaHF(:,Nt*0.75),'LineWidth',3);
plot(xi,etaHF(:,Nt*1),'-o','LineWidth',3);
legend('\tau=0.25\tau_{total}', '\tau=0.5\tau_{total}', '\tau=0.75\tau_{total}','\tau=\tau_{total}')
xlabel('\xi'); ylabel('{\eta}_{HF}');
prop_plots

figure
plot(tau,etaHF(Nx*0.25,:),'--','LineWidth',3); hold on
plot(tau,etaHF(Nx*0.5,:),'LineWidth',3); plot(tau,etaHF(Nx*0.75,:),'LineWidth',3);
plot(tau,etaHF(Nx*1,:),'-o','LineWidth',3);
legend('\xi=0.25\xi_{total}', '\xi=0.5\xi_{total}', '\xi=0.75\xi_{total}','\xi=\xi_{total}')
xlabel('\tau'); ylabel('{\eta}_{HF}');
prop_plots


%%========================================================
%% QoI: Vcell
etaHF_xi0 = etaHF(1,:);
etaHF_xi1 = etaHF(end,:);
etaLF_xi0 = etaLF(1,:);
etaLF_xi1 = etaLF(end,:);

V_elecHF = ((1+2*gamma)/(1+gamma))*etaHF_xi1 - (gamma/(1+gamma))*etaHF_xi0 - (gamma/(1+gamma)^2)*I;
V_elecLF = ((1+2*gamma)/(1+gamma))*etaLF_xi1 - (gamma/(1+gamma))*etaLF_xi0 - (gamma/(1+gamma)^2)*I;

B = (Ls/kappa_s)*(kappa*sigma/(L*(kappa+sigma)));

V_cellHF = 1 - 0.5*B*I - V_elecHF;
V_cellLF = 1 - 0.5*B*I - V_elecLF;

figure
plot(tau,V_cellHF,'--b','LineWidth',3); hold on
plot(tau,V_cellLF,'r','LineWidth',3);
legend('HF model', 'LF model')
xlabel('\tau'); ylabel('V_{cell}');
prop_plots

figure
plot(I,V_elecHF,'--b','LineWidth',3); hold on
plot(I,V_elecLF,'r','LineWidth',3);
legend('HF model', 'LF model')
xlabel('I'); ylabel('V^*_{cell}');
prop_plots

%% Error in Velectrod

error_Velec = V_elecHF - V_elecLF;

figure
plot(tau,error_Velec,'b','LineWidth',3); hold on
xlabel('\tau'); ylabel('\epsilon_{Velec}');
prop_plots

%%========================================================
%% Save in data file
M = [tau' I' error_Velec'];
csvwrite('error_in_qoi.txt',M);


