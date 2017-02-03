%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve the HF model evolution equation (Parabolic PDE)
%%%
%%%    d eta/ d tau = d2eta/dxi2 + rho(xi,tau)
%%%    BCs:     deta/dxi = alpha  at xi=0
%%%             deta/dxi = beta   at xi=1
%%%    ICs:     eta = 0    at tau=0
%%%
%%% Crank Nicolson time-stepping, and Finite Difference in space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = fdm(t_0,t_f,x_0,x_f, Nt,Nx, alpha, beta, eta0)

% t_0   is the initial time,
% t_f   is the final time,
% Nx    is the number of mesh-points,
% Nt    is the number of time steps.
% rho   is the source terms (Nx*Nt materix)
% alpha is the BC at xi=0 (1*Nt vector)
% beta  is the BC at xi=1 (1*Nt vector)


%% Define the mesh in Space and Time
dx = (x_f - x_0)/(Nx-1);
x = 0:dx:(x_f - x_0);
x = x';
dt = (t_f-t_0)/(Nt-1);
t = t_0:dt:t_f;

% define the ratio r
r = 0.5*dt/dx^2;

u = zeros(Nx,Nt);

%% Impose IC
% for i=1:Nx
% %     u(i,1) = 1-cos(2*pi*i/Nx);
     u(:,1) = eta0;
% end

% construct A of spatial system
Ax = zeros(Nx,Nx);
Ax(1,1) = -2;
Ax(1,2) = 2;
for i=2:Nx-1
    Ax(i,i-1) = 1;
    Ax(i,i)   = -2;
    Ax(i,i+1) = 1;
end
Ax(Nx,Nx)   = -2;
Ax(Nx,Nx-1) = 2;

Ax = (1/dx^2)*Ax;


% construct A and B of full system
% AU(n+1) = BU(n) + C 

A = eye(Nx,Nx) - 0.5*dt*Ax;

B = eye(Nx,Nx) + 0.5*dt*Ax;


for j=1:Nt-1
    C(1:Nx) = 0;    
    C(1) = -0.5*dt*(2*alpha(j)/dx+2*alpha(j+1)/dx);
    C(Nx) = -0.5*dt*(-2*beta(j)/dx-2*beta(j+1)/dx);
    
    BUmatrix = B*u(:,j);
    for i=1:Nx
        RHS(i) = BUmatrix(i)   + C(i);
    end    
    u(:,j+1) = A\RHS';
end