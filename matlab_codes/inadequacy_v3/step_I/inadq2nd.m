
function sse = inadq1st(x,tau,eps_exact, I, Nt)




lambda = x(1);
mu = x(2);
alpha = x(3);
beta = x(4);
rho = x(5);
eps1_zero = x(6);
eps2_zero = x(7);

dtau = tau(2) - tau(1);
Epsilon1(1) = eps1_zero;
Epsilon2(1) = eps2_zero;

for i = 1 : Nt-2
    Epsilon1(i+1) = Epsilon1(i) + Epsilon2(i)*dtau;
    Epsilon2(i+1) = Epsilon2(i) + ...
        dtau*(alpha*I(i) - mu^2*Epsilon1(i) - lambda^2*Epsilon2(i)) + ...
        beta*(I(i+1)-I(i)) + (rho/dtau)*(I(i+2)-2*I(i+1)+I(i));
end


eps_exact_shor = eps_exact(1:Nt-1);

sse = sum((eps_exact_shor - Epsilon1).^2);
