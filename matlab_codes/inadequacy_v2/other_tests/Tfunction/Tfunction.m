%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script to plot inadequecy representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
set(0,'DefaultFigureWindowStyle','docked')


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
% Ins = Iamp*square(6*pi*tau');
Ins = Iamp*sin(6*pi*tau');
% Ins(181:359) = -Ins(181:359);
I = Ins';



%% T(tau) evolution
dtau = tau(2) - tau(1);
T(1:Ntau) = 0;


for i = 1 : Ntau-1
    
    T(i+1) = T(i) + (tau(i) - T(i))*dtau*abs((I(i+1)-I(i))/dtau);
    
end

%%


figure
plot(tau,T,'r','LineWidth',3); hold on
xlabel('\tau'); ylabel('T(\tau)');
prop_plots
axis square

figure
plot(tau,I,'b','LineWidth',3);
xlabel('\tau'); ylabel('I(\tau)');
prop_plots
axis square
