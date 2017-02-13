%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script to plot inadequecy representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
set(0,'DefaultFigureWindowStyle','docked')


%% Loading HF data
A = csvread('error_in_qoi.txt');
tau       = A(:,1);
I         = A(:,2);
eps_exact = A(:,3); 

% figure
% plot(tau,eps_exact,'LineWidth',3)
% xlabel('\tau'); ylabel('\epsilon');
% prop_plots


%% Inadequacy Parameters
Nsamp = 20;
lambda_mean = 15;
c = 15;
alpha = 0.32;
beta = 15;

% initial conditions
epsODE(1,1:Nsamp) = 0.02*randn(1,Nsamp)+0.14;
lambda(1,1:Nsamp) = 20.0;


%% Inadequecy evolution
Ntau = size(tau,1);
dtau = tau(2) - tau(1);

% Wiener process
dW = sqrt(dtau)*randn(Nsamp,Ntau);   % increments
dW = dW';

for i = 1 : Ntau-1
    for j = 1 : Nsamp
    
    lambda(i+1,j) = lambda(i,j) - c*(lambda(i,j)-lambda_mean)*dtau + beta*dW(i,j);
    epsODE(i+1,j) = epsODE(i,j) - lambda(i,j)*epsODE(i,j)*dtau + alpha*(I(i+1)-I(i));

    end
end

%%
figure
plot(tau,lambda,'LineWidth',1)
xlabel('\tau'); ylabel('\lambda');
prop_plots

figure
plot(tau,eps_exact,'b','LineWidth',3); hold on
plot(tau,-epsODE,'--','LineWidth',2)
xlabel('\tau'); ylabel('Epsilon');
prop_plots


%% confidence plot
Meps = mean(-epsODE');
Leps = min(-epsODE');
Ueps = max(-epsODE');

figure
p = plot(tau,Leps,tau, Ueps,tau,Meps);
YLIM = get(gca,'YLim'); delete(p);
a1 = area(tau,Ueps,min(YLIM)); 
hold on
set(a1,'LineStyle','none');     set(a1,'FaceColor',[0.9 0.9 0.9]);
a2 = area(tau,Leps,min(YLIM)); 
set(a2,'LineStyle','none');     set(a2,'FaceColor',[1 1 1]);
p_mean = plot(tau,Meps,'--b','LineWidth',2);
p_exact = plot(tau, eps_exact, 'r','LineWidth',2); hold on
legend([p_exact p_mean a1],{'exact','mean','95% confidence boundaries'})
xlabel('Normalized Time \tau'); ylabel('Error in QoI');
prop_plots
daspect([1 2 1])
