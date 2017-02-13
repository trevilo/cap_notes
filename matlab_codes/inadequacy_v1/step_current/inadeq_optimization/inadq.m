
function sse = inadq(x,tau,eps_exact, I, Nt)


c = x(1);
lambda_mean = x(2);
alpha = x(3);

dtau = tau(2) - tau(1);
lambda(1) = 20.0;
Epsilon(1) = 0.12;

for i = 1 : Nt-1
    
    lambda(i+1) = lambda(i) - c*(lambda(i)-lambda_mean)*dtau;
    Epsilon(i+1) = Epsilon(i) - (lambda(i)*Epsilon(i)*dtau) + alpha*(I(i+1)-I(i));

end

sse = sum((-eps_exact - Epsilon).^2);
