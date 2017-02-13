%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script for calibrating inadequecy ODE
%%% against HF data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
set(0,'DefaultFigureWindowStyle','docked')


%% Loading HF data
A = csvread('../hf_lf_models/error_in_qoi.txt');
tau_org       = A(:,1);
I_org         = A(:,2);
eps_exact_org = A(:,3);

tau = tau_org(181:end) - tau_org(181);
I = I_org(181:end);
eps_exact = eps_exact_org(181:end);

plot(eps_exact)

% eps_exact = csvread('error_cycl.txt');
% 
% Nt = 800;  %--- number of time steps
% tau = linspace(0,0.743319442711514 , Nt);  %--- sec
% Ins = 0.410044212529724*square(6*pi*tau');
% I = Ins';
% plot(I)

%% optmization (SQ function)
funSQ = @(x)inadqSQ(x,tau,eps_exact, I);

x0 = 0.2;
bestxSQ = fminsearch(funSQ,x0)


lambda_mean_b = bestxSQ;
dtau = tau(2) - tau(1);
Epsilon0 = 0.2397;

EpsSQfit = Epsilon0 - 2*lambda_mean_b*sqrt(tau);


%% optmization (EXP function)
funEXP = @(xEXP)inadqEXP(xEXP,tau,eps_exact, I);

x0EXP = 2.5;
bestxEXP = fminsearch(funEXP,x0EXP)


lambda_mean_bEXP = bestxEXP;
dtau = tau(2) - tau(1);
Epsilon0 = 0.2397;

EpsEXPfit = Epsilon0*exp(-2*lambda_mean_bEXP*sqrt(tau));


%%
plot(tau+0.1851,eps_exact,'b','LineWidth',4);
hold on
plot(tau+0.1851,EpsSQfit,'--r','LineWidth',4);
plot(tau+0.1851,EpsEXPfit,':k','LineWidth',4);
xlabel('\tau'); ylabel('\epsilon');
legend('HF data','fitted option I','fitted option II')
axis square
prop_plots