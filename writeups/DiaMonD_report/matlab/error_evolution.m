%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve the error evolution equation (Parabolic PDE)
%%% 
%%%    d eps/ d tau = d2eps/dxi2 - rho(xi,tau)
%%%    BCs:     deps/dxi= 0  at xi=0,1
%%%    ICs:     GIVEN    at tau=0
%%%
%%% Crank Nicolson time-stepping, and Finite Difference in space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = error_evolution(t_0,t_f,x_0,x_f, Nt,Nx, rho, ICond, alpha, beta)

% this solves the equation u_t = u_xx - F with 
%       Neumann boundary conditions
%       Initial conditions : u(t=0) = ICond 
% t_0 is the initial time, t_f is the final time, 
% Nx is the number of mesh-points, and Nt is the number of time steps. 


% define the mesh in space
dx = (x_f - x_0)/(Nx-1);
x = 0:dx:(x_f - x_0);
x = x';

% define the mesh in time
dt = (t_f-t_0)/(Nt-1);
t = t_0:dt:t_f;

% define the ratio r
r = dt/dx^2;

% Impose IC
for i=1:Nx
  u(i,1) = ICond(i,1);
end 

%% Creating materixes [A]{u_new} = [B]{u_old} + {C_old}
% for internal points, have
%    u_new(j) = u_old(j) + r/2*(u_old(j+1)-2*u_old(j)+u_old(j-1))
% for the two end-points, have
%    u_new(1) = u_old(1) + r/2*(2 u_old(2)-2*u_old(1))
%    u_new(N) = u_old(N) + r/2*( -2*u_old(N)+2 u_old(N-1))     
% clearly the endpoints are redundant: u(1)= u(N) at all times.  I just
% kept them around for plotting convenience.  

% define the matrix A which has to be inverted at every time-step.
%   u_new(1) - u_old(1) = dt/h^2 (2 u_new(2)-2*u_new(1))
%   u_new(i) - u_old(i) = dt/h^2 (u_new(i+1)-2*u_new(i)+u_new(i-1)
%   u_new(N) - u_old(N) = dt/h^2 ( -2*u_new(N)+ 2 u_new(N-1))


A = zeros(Nx-1,Nx-1);
A(1,1) = 1+r;
A(1,2) = -r;
for i=2:Nx-1
  A(i,i-1) = -r/2;
  A(i,i) = 1+r;
  A(i,i+1) = -r/2;
end
A(Nx,Nx) = 1+r;
A(Nx,Nx-1) = -r;

% for j=1:Nt
%   RHS(1) = u(1,j) + r/2*(2*u(2,j)-2*u(1,j));% - rho(1,j)*dt;
%   for i=2:Nx-1
%     RHS(i) = u(i,j) + r/2*(u(i+1,j)-2*u(i,j)+u(i-1,j)) - rho(i,j)*dt;
%   end
%   RHS(Nx) = u(Nx,j) + r/2*(-2*u(Nx,j)+2*u(Nx-1,j));% - rho(Nx,j)*dt;
%   u(:,j+1) = tridiag(A,RHS);  
% end

for j=1:Nt-1
  RHS(1) = u(1,j) + r/2*(2*u(2,j)-2*u(1,j)) + rho(1,j)*dt  - r*alpha(j)*dx;
  for i=2:Nx-1
    RHS(i) = u(i,j) + r/2*(u(i+1,j)-2*u(i,j)+u(i-1,j)) + rho(i,j)*dt;
  end
  RHS(Nx) = u(Nx,j) + r/2*(-2*u(Nx,j)+2*u(Nx-1,j)) + rho(Nx,j)*dt + r*beta(j)*dx;
%    u(:,j+1) = tridiag(A,RHS);
   u(:,j+1) = inv(A)*RHS';
end

% 
% surf(t,x,u)