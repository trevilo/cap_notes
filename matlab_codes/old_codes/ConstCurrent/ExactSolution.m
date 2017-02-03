%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Exact Solution of eta
%%%     from Srinivasan and Wedner 1999
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all;

gamma = 0;
I = 1;
xi = 0:1/100:1;
tau = [0.65; 0.4; 0.2; 0.1; 0.04; 0.02];
syms n

for j = 1:size(tau,1)
for i = 1: size(xi,2)
    S = symsum((((-1)^n+gamma)/n^2)*cos(n*pi*xi(i))*exp(-(n)^2*pi^2*tau(j)), n, 1, 100);
    eta(i,j) = I*tau(j) + I*((3*xi(i)^2-1)/(6*(1+gamma))) + I*(gamma*(3*xi(i)^2-6*xi(i)+2)/(6*(1+gamma))) - ...
        2*I/(pi^2*(1+gamma)) * S;
end
end

plot(xi,eta)
