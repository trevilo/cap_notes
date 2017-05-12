
function sse = inadq1st(x,tau,eps_exact, I, Nt)


lambda = x(1);
beta = x(2);
eps_zero = x(3);


dtau = tau(2) - tau(1);

Epsilon(1) = eps_zero;

for i = 1 : Nt-1
    Epsilon(i+1) = Epsilon(i) - (lambda^2*Epsilon(i)*dtau) + beta*dtau*I(i);
end

sse = sum((eps_exact - Epsilon).^2);
