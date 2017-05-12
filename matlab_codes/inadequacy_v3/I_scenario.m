%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Script to generate HF data
%%%          for constant current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all


%% Parameters
kappa = 0.0195174;   %--- sec/m
sigma = 52.1;   %--- sec/m
gamma = kappa/sigma;
L = 50e-6;  %--- m
C = 0.03134;  %--- F/m2
a = 4.19956e7/C;   %--- m
V0 = 1.25;   %--- volt
Iunscaled = 200;   %--- Amp/m^2
Ls = 25e-6;  %--- m
kappa_s = 0.0311627;   %--- sec/m

Nt = 5000;  %--- number of time steps
Nx = 700;  %--- number of space steps 

%% Domain
t = linspace(0.001,5 , Nt);  %--- sec
x = linspace(0,L , Nx);

% convert to nondimesional values
tau = (kappa*sigma/(kappa+sigma))*(1/(a*C*L^2)) .* t;
xi = x ./ L;
Nxi = Nx ; Ntau = Nt;

%% Applied Current
t_conv = a*C*L^2*(kappa+sigma)/(kappa*sigma);

I_const = 0.5* (Iunscaled*(L/V0)*((kappa+sigma)/(sigma*kappa))) * ones(1,Nt);
I_const(1) = 0;

Iamp = (Iunscaled*(L/V0)*((kappa+sigma)/(sigma*kappa)));
Insstep = Iamp*square(6*pi*tau');
I = Insstep';
I_step = smooth(Insstep,'moving');
I_step(1) = 0;

Ins_sin = 1.5*Iamp*sin(5*pi*tau');
I_sin = Ins_sin';
%I = smooth(Ins_sin,'moving');

%%
figure
% plot(tau,I_const,'-b','LineWidth',3); hold on
plot(tau,I_step,'-g','LineWidth',3); hold on
plot(tau,I_sin,'-b','LineWidth',3); 
legend( 'step', 'sinosoidal')
xlabel('$$\tau$$','Interpreter','LaTex'); ylabel('$I$','Interpreter','LaTex');
prop_plots
