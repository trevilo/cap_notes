%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script to provide data for calibrating 
%%%         model inadequecy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
set(0,'DefaultFigureWindowStyle','docked')

error_exact_Velec = csvread('exact_error.txt');

%% Inadequacy Parameters
lambda_mean = 12;
c = 12;
beta = 30;

Epsilon(1) = 0.14;
lambda(1) = 20.0;

kappa = 0.0195174;   %--- sec/m
sigma = 52.1;   %--- sec/m
gamma = kappa/sigma;
L = 50e-6;  %--- m
C = 0.03134;  %--- F/m2
a = 4.19956e7/C;   %--- m

%% Domain
L  = 50e-6;
Nt = 700;  %--- number of time steps
Nx = 700;  %--- number of space steps 

t = linspace(0,3 , Nt);  %--- sec
x = linspace(0,L , Nx);

% convert to nondimesional values
tau = (kappa*sigma/(kappa+sigma))*(1/(a*C*L^2)) .* t;
xi = x ./ L;
Nxi = Nx ; Ntau = Nt;
t_conv = a*C*L^2*(kappa+sigma)/(kappa*sigma);

dtau = tau(2) - tau(1);

%% Inadequecy evolution
% randn('state',100)
dW = sqrt(dtau)*randn(1,Ntau);   % increments

for i = 1 : Ntau-1
    
    lambda(i+1) = lambda(i) - c*(lambda(i)-lambda_mean)*dtau + beta*dW(i);
    Epsilon(i+1) = Epsilon(i) - lambda(i)*Epsilon(i)*dtau;

end

figure
plot(tau,lambda,'--r','LineWidth',3)
xlabel('\tau'); ylabel('\lambda');
prop_plots

figure
plot(tau,abs(error_exact_Velec),'b','LineWidth',3); hold on
plot(tau,Epsilon,'-k','LineWidth',2)
xlabel('\tau'); ylabel('Epsilon');
prop_plots

