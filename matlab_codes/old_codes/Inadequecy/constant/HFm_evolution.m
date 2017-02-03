%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve the HF model evolution equation (Parabolic PDE)
%%% 
%%%    d eta/ d tau = d2eta/dxi2 + rho(xi,tau)
%%%    BCs:     deta/dxi = alpha  at xi=0
%%%             deta/dxi = beta   at xi=1
%%%    ICs:     eps = 0    at tau=0
%%%
%%% Crank Nicolson time-stepping, and Finite Difference in space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = HFm_evolution(t_0,t_f,x_0,x_f, Nt,Nx, rho, alpha, beta)

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
for i=1:Nx
%   u(i,1) = sin(pi*x(i)/(x_f - x_0));
  u(i,1) = 0;  
end 

%% Neumann BCs
% % at xi = 0
% alpha = 10*t;
% 
% % at xi = 1
% beta = 0*t;

%% Creating materixes [A]{u_new} = [B]{u_old} + {C_old}
% for internal points, have
%    u_new(j) = u_old(j) + r/2*(u_old(j+1)-2*u_old(j)+u_old(j-1))
% for the two end-points, have
%    u_new(1) = u_old(1) + r/2*(2 u_old(2)-2*u_old(1))
%    u_new(N) = u_old(N) + r/2*( -2*u_old(N)+2 u_old(N-1))     

% define the matrix A which has to be inverted at every time-step.
%   u_new(1) - u_old(1) = dt/h^2 (2 u_new(2)-2*u_new(1))
%   u_new(i) - u_old(i) = dt/h^2 (u_new(i+1)-2*u_new(i)+u_new(i-1)
%   u_new(N) - u_old(N) = dt/h^2 ( -2*u_new(N)+ 2 u_new(N-1))


A = zeros(Nx-1,Nx-1);
A(1,1) = 1+r;
A(1,2) = -r;
for i=2:Nx-1
  A(i,i-1) = -r;
  A(i,i) = 1+2*r;
  A(i,i+1) = -r;
end
A(Nx,Nx) = 1+r;
A(Nx,Nx-1) = -r;

B = zeros(Nx-1,Nx-1);
B(1,1) = 1-r;
B(1,2) = r;
for i=2:Nx-1
  B(i,i-1) = r;
  B(i,i) = 1-2*r;
  B(i,i+1) = r;
end
B(Nx,Nx) = 1-r;
B(Nx,Nx-1) = r;


for j=1:Nt-1
   BUmatrix = B*u(:,j); 
   RHS(1) = BUmatrix(1)   + rho(1,j)*dt  - 2*r*alpha(j)*dx;
   for i=2:Nx-1
       RHS(i) = BUmatrix(i)   + rho(i,j)*dt;
   end
   RHS(Nx) = BUmatrix(Nx) + rho(Nx,j)*dt + 2*r*beta(j)*dx;
   u(:,j+1) = tridiag(A,RHS);
end

% for j=1:Nt-1
%   RHS(1) = u(1,j) + r/2*(2*u(2,j)-2*u(1,j)) + rho(1,j)*dt - r*alpha(j)*dx;
%   for i=2:Nx-1
%     RHS(i) = u(i,j) + r/2*(u(i+1,j)-2*u(i,j)+u(i-1,j)) + rho(i,j)*dt;
%   end
%   RHS(Nx) = u(Nx,j) + r/2*(-2*u(Nx,j)+2*u(Nx-1,j)) + rho(Nx,j)*dt + r*beta(j)*dx ;
%   u(:,j+1) = tridiag(A,RHS);  
% end

% surf(t,x,u)
