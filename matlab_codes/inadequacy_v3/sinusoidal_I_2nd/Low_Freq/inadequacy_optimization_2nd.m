%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script for calibrating inadequecy ODE
%%% against HF data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc


%% Loading HF data
A = csvread('error_in_qoi_sinI.txt');
tau       = A(:,1);
I         = A(:,2);
eps_exact = A(:,3); 

Nt = size(tau,1);

%% optmization
fun = @(x)inadq2nd(x,tau',eps_exact',I, Nt);
 
x0 = [8.8653   26.6089   36.3662  -15.3801   -0.1414   -1.9198];
bestx = fminsearch(fun,x0)


lambda_b = bestx(1);
mu_b = bestx(2);
alpha_b = bestx(3);
beta_b = bestx(4);
rho_b = bestx(5);
eps1_zero_b = 0;
eps2_zero_b = bestx(6);

dtau = tau(2) - tau(1);
Epsilonfit1(1) = eps1_zero_b;
Epsilonfit2(1) = eps2_zero_b;

for i = 1 : Nt-2
    Epsilonfit1(i+1) = Epsilonfit1(i) + Epsilonfit2(i)*dtau;
    Epsilonfit2(i+1) = Epsilonfit2(i) + ...
        dtau*(alpha_b*I(i) - mu_b^2*Epsilonfit1(i) - lambda_b^2*Epsilonfit2(i)) + ...
        beta_b*(I(i+1)-I(i)) + (rho_b/dtau)*(I(i+2)-2*I(i+1)+I(i));
end


plot(tau,eps_exact,'b','LineWidth',3);
hold on
plot(tau(1:Nt-1),Epsilonfit1,'--g','LineWidth',3);
xlabel('\tau'); ylabel('\epsilon');
legend('HF data','fitted 2nd ODE')
axis square
prop_plots