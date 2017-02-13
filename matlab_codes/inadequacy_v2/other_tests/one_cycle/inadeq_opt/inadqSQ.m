
function sse = inadqSQ(x,tau,eps_exact, I)


lambda_mean = x;


dtau = tau(2) - tau(1);
Epsilon0 = 0.2397;

EpsSQ = Epsilon0 - 2*lambda_mean*sqrt(tau);


sse = sum((eps_exact - EpsSQ).^2);
