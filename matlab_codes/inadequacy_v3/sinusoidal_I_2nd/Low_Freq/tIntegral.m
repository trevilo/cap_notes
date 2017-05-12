function V = tIntegral(tau)
kappa = 0.0195174;   %--- sec/m
sigma = 52.1;   %--- sec/m

L = 50e-6;  %--- m
C = 0.03134;  %--- F/m2

V0 = 1.25;   %--- volt
Iunscaled = 200;   %--- Amp/m^2

Iamp = (Iunscaled*(L/V0)*((kappa+sigma)/(sigma*kappa)));
% Ins = Iamp*square(6*pi*tau');
fun = @(t,tau) Iamp*sin(5*pi*tau').*k(t,tau);
V = integral(@(t)fun(t,tau),0,tau);
end