%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script for calibrating inadequecy ODE
%%% against HF data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc


%% Loading HF data
A = csvread('error_in_qoi_stepI.txt');
tau       = A(:,1);
I         = A(:,2);
eps_exact = A(:,3); 

Nt = size(tau,1);

%% optmization
fun = @(x)inadq1st(x,tau',eps_exact',I, Nt);

x0 = [4.2    0.8    0.13    0.3];
bestx = fminsearch(fun,x0)


lambda_b = bestx(1);
beta_b = bestx(2);
eps_zero_b = bestx(3); 
alpha_b = bestx(4); 

dtau = tau(2) - tau(1);
Epsilonfit(1) = eps_zero_b;

for i = 1 : Nt-1
    Epsilonfit(i+1) = Epsilonfit(i) - (lambda_b^2*Epsilonfit(i)*dtau) + beta_b*dtau*I(i) + alpha_b*(I(i+1)-I(i));
end


plot(tau,eps_exact,'b','LineWidth',4);
hold on
plot(tau,Epsilonfit,'--r','LineWidth',3);
xlabel('\tau'); ylabel('\epsilon');
legend('HF data','fitted 1st order ODE')
axis square
prop_plots
