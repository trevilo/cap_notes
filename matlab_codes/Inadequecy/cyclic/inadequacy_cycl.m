%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script to provide data for calibrating 
%%%         model inadequecy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
set(0,'DefaultFigureWindowStyle','docked')

error_exact_Velec = csvread('error_cycl.txt');
V_elecHF = csvread('V_elecHF.txt');
V_elecLF = csvread('V_elecLF.txt');

%% Inadequacy Parameters
Nsamp = 17;
lambda_mean = 15;
c = 15;
alpha = 0.32;
beta = 15;

% Epsilon(1,1:Nsamp) = 0.14;
Epsilon(1,1:Nsamp) = 0.02*randn(1,Nsamp)+0.14;
lambda(1,1:Nsamp) = 20.0;

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
Time = 4; %--- sec

Nt = 800;  %--- number of time steps
Nx = 500;  %--- number of space steps 


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

% I = (Iunscaled*(L/V0)*((kappa+sigma)/(sigma*kappa))) * ones(1,Nt);


dtau = tau(2) - tau(1);

%% Inadequecy evolution
% randn('state',100)
dW = sqrt(dtau)*randn(Nsamp,Ntau);   % increments
dW = dW';

for i = 1 : Ntau-1
    for j = 1 : Nsamp
    
    lambda(i+1,j) = lambda(i,j) - c*(lambda(i,j)-lambda_mean)*dtau + beta*dW(i,j);
    Epsilon(i+1,j) = Epsilon(i,j) - lambda(i,j)*Epsilon(i,j)*dtau + alpha*(I(i+1)-I(i));

    end
end

%%
% figure
% plot(tau,lambda,'--r','LineWidth',3)
% xlabel('\tau'); ylabel('\lambda');
% prop_plots
% 
% figure
% plot(tau,error_exact_Velec,'b','LineWidth',3); hold on
% plot(tau,-Epsilon,'--','LineWidth',2)
% xlabel('\tau'); ylabel('Epsilon');
% prop_plots
%%
figure
plot(tau,V_elecHF,'r','LineWidth',3); hold on
plot(tau,V_elecLF,'--b','LineWidth',3);
xlabel('Normalized Time \tau'); ylabel('QoI: V^{cell}');
legend ('V_{HF}','V_{LF}')
prop_plots
daspect([1 2 1])

%% confidence plot
Meps = mean(-Epsilon');
Leps = min(-Epsilon');
Ueps = max(-Epsilon');

figure
p = plot(tau,Leps,tau, Ueps,tau,Meps);
YLIM = get(gca,'YLim'); delete(p);
a1 = area(tau,Ueps,min(YLIM)); 
hold on
set(a1,'LineStyle','none');     set(a1,'FaceColor',[0.9 0.9 0.9]);
a2 = area(tau,Leps,min(YLIM)); 
set(a2,'LineStyle','none');     set(a2,'FaceColor',[1 1 1]);
p_mean = plot(tau,Meps,'--b','LineWidth',2);
p_exact = plot(tau,error_exact_Velec,'r','LineWidth',2); hold on
legend([p_exact p_mean a1],{'exact','mean','95% confidence boundaries'})
xlabel('Normalized Time \tau'); ylabel('Error in QoI');
prop_plots
daspect([1 2 1])

%% plot QoI correction
for i = 1 : Nsamp
    Vlf_corrected(i,:) = V_elecLF(:)-Epsilon(:,i);
end
Mv = mean(Vlf_corrected);
Lv = min(Vlf_corrected);
Uv = max(Vlf_corrected);

figure
p = plot(tau,Lv,tau, Uv,tau,Mv);
YLIM = get(gca,'YLim'); delete(p);
a1 = area(tau,Uv,min(YLIM)); 
hold on
set(a1,'LineStyle','none');     set(a1,'FaceColor',[0.9 0.9 0.9]);
a2 = area(tau,Lv,min(YLIM)); 
set(a2,'LineStyle','none');     set(a2,'FaceColor',[1 1 1]);
p_mean = plot(tau,Mv,'--b','LineWidth',2);
p_exact = plot(tau,V_elecHF,'r','LineWidth',3);
xlabel('Normalized Time \tau'); ylabel('QoI: V^{cell}');
legend([p_exact p_mean a1],{'V_{HF}','V_{LF}+error','95% confidence boundaries'})
prop_plots
daspect([1 2 1])
