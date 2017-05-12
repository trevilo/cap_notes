%% Kernel term of the solution

function integrand = k(t,tau)
kappa = 0.0195174;   %--- sec/m
sigma = 52.1;   %--- sec/m
gamma = kappa/sigma;

integrand = 1;
for n=1:50000
    integrand = integrand + 2*(2*gamma/(1+gamma)^2*(-1)^n+(1+gamma^2)/(1+gamma)^2)...
        *exp(-n^2*pi^2*(tau-t));
end
end
