clear all
close all
clc


%% Loading HF data
A = csvread('error_in_qoi.txt');
tau       = A(:,1);
I         = A(:,2);
eps_exact = A(:,3); 

Nt = size(tau,1);

%% inadequacy

%-- 1 term
C1 = 0.1033;
t1 = 0.0855;
epsI = C1*exp(-tau./t1);

%-- 2 terms
C1 = 0.0362;
t1 = 0.0132;
C2 = 0.087;
t2 = 0.0996;
 
epsII = C1*exp(-tau./t1) + C2*exp(-tau./t2);


%% Plot
figure
plot(tau,eps_exact,':k','LineWidth',5);
hold on
plot(tau,epsI,'b','LineWidth',3);
plot(tau,epsII,'g','LineWidth',3);
xlabel('\tau'); ylabel('\epsilon');
legend('HF data','fitted ODE (v1)','fitted ODE (v2)')
axis square
prop_plots