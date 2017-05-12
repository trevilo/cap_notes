
function sse = inadq1st(x,tau,eps_exact, I, Nt)




lambda = x(1);
mu = x(2);
alpha = x(3);
eps1_zero = x(4);
eps2_zero = x(5);

dtau = tau(2) - tau(1);
Epsilon1(1) = eps1_zero;
Epsilon2(1) = eps2_zero;

for i = 1 : Nt-1
    Epsilon1(i+1) = Epsilon1(i) + Epsilon2(i)*dtau;
    Epsilon2(i+1) = Epsilon2(i) + ...
        dtau*(alpha*I(i) - mu^2*Epsilon1(i) - lambda^2*Epsilon2(i));
end

sse = sum((eps_exact - Epsilon1).^2);
