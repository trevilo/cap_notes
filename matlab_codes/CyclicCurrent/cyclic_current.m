%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script to solve the simple example
%%%     for model inadequecy
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
Time = 4; %--- sec

Nt = 800;  %--- number of time steps
Nx = 800;  %--- number of space steps 


%% Domain
t = linspace(0,Time , Nt);  %--- sec
x = linspace(0,L , Nx);

% convert to nondimesional values
tau = (kappa*sigma/(kappa+sigma))*(1/(a*C*L^2)) .* t;
xi = x ./ L;
Nxi = Nx ; Ntau = Nt;

%% Applied Current
t_conv = a*C*L^2*(kappa+sigma)/(kappa*sigma);

% I = (Iunscaled*(L/V0)*((kappa+sigma)/(sigma*kappa))) * ones(1,Nt);

%== Step Function
I = (Iunscaled*(L/V0)*((kappa+sigma)/(sigma*kappa))) * ones(1,Nt);

for i = 1*Nt/4 : 2*Nt/4
    I(i) = -I(i);
end
for i = 3*Nt/4 : 4*Nt/4
    I(i) = -I(i);
end

figure(1)
plot(tau,I,'LineWidth',3);
axis square
set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'String','\tau','FontSize',32,'FontWeight','bold','FontName','Times')
set(get(gca,'YLabel'),'String','I','FontSize',32,'FontWeight','bold','FontName','Times')
set(gcf,'Position',[1 1 round(1000) round(1000)])
set(get(gca,'Title'),'String','Cyclic Current','FontSize',32,'FontWeight','bold','FontName','Times')


%%========================================================
%% Solution of LF model
eta_bar = tau.*I;

for i = 1 : size(tau,2) % time
    for j = 1 : size(xi,2) % space
        etaLF(j,i) = (I(i)/2)*xi(j)^2 - ((I(i)*gamma)/(1+gamma))*xi(j) + eta_bar(i) - I(i)/6 + (I(i)*gamma)/(2*I(i)+2*gamma);
    end
end

figure
surf(tau, xi, etaLF,'EdgeColor','none')
% axis square
set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'String','\tau','FontSize',42,'FontWeight','bold','FontName','Times')
set(get(gca,'YLabel'),'String','\xi','FontSize',42,'FontWeight','bold','FontName','Times')
set(get(gca,'ZLabel'),'String','{\eta}_{LF}','FontSize',42,'FontWeight','bold','FontName','Times')
set(get(gca,'Title'),'String','Cyclic Current','FontSize',32,'FontWeight','bold','FontName','Times')
set(gcf,'Position',[1 1 round(1000) round(1000)])

figure
plot(xi,etaLF(:,Nt*0.25),'--','LineWidth',3); hold on
plot(xi,etaLF(:,Nt*0.5),'LineWidth',3);
plot(xi,etaLF(:,Nt*0.75),'LineWidth',3);
plot(xi,etaLF(:,Nt*1),'-o','LineWidth',3);
legend('\tau=0.25\tau_{total}', '\tau=0.5\tau_{total}', '\tau=0.75\tau_{total}','\tau=\tau_{total}')
axis square
set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'String','\xi','FontSize',42,'FontWeight','bold','FontName','Times')
set(get(gca,'YLabel'),'String','{\eta}_{LF}','FontSize',42,'FontWeight','bold','FontName','Times')
set(get(gca,'Title'),'String','Cyclic Current','FontSize',32,'FontWeight','bold','FontName','Times')
set(gcf,'Position',[1 1 round(1000) round(1000)])

%%========================================================
%% Solution of HF model
% 
% % computing Neumann BCs:
%-- alpha at xi=0
alpha = -I.*(gamma/(1+gamma));
%-- beta at xi=1
beta  =  I.*(1/(1+gamma));

% solving HF pde+bc using finite difference
etaHF = HFm_evolution(tau(1),tau(end),xi(1),xi(end), Ntau,Nxi, zeros(Nxi,Ntau), alpha, beta);

figure
surf(tau, xi, etaHF,'EdgeColor','none')
% axis square
set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'String','\tau','FontSize',32,'FontWeight','bold','FontName','Times')
set(get(gca,'YLabel'),'String','\xi','FontSize',32,'FontWeight','bold','FontName','Times')
set(get(gca,'ZLabel'),'String','{\eta}_{HF}','FontSize',32,'FontWeight','bold','FontName','Times')
set(get(gca,'Title'),'String','Cyclic Current','FontSize',32,'FontWeight','bold','FontName','Times')
set(gcf,'Position',[1 1 round(1000) round(1000)])

figure
plot(xi,etaHF(:,Nt*0.25),'--','LineWidth',3); hold on
plot(xi,etaHF(:,Nt*0.5),'LineWidth',3);
plot(xi,etaHF(:,Nt*0.75),'LineWidth',3);
plot(xi,etaHF(:,Nt*1),'-o','LineWidth',3);
legend('\tau=0.25\tau_{total}', '\tau=0.5\tau_{total}', '\tau=0.75\tau_{total}','\tau=\tau_{total}')
axis square
set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'String','\xi','FontSize',42,'FontWeight','bold','FontName','Times')
set(get(gca,'YLabel'),'String','{\eta}_{HF}','FontSize',42,'FontWeight','bold','FontName','Times')
set(get(gca,'Title'),'String','Cyclic Current','FontSize',32,'FontWeight','bold','FontName','Times')
set(gcf,'Position',[1 1 round(1000) round(1000)])

%%========================================================
%% Error evolution = HF model - LF model

eps_exact = etaHF - etaLF;

figure
surf(tau, xi, eps_exact,'EdgeColor','none')
% axis square
set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'String','\tau','FontSize',32,'FontWeight','bold','FontName','Times')
set(get(gca,'YLabel'),'String','\xi','FontSize',32,'FontWeight','bold','FontName','Times')
set(get(gca,'ZLabel'),'String','\epsilon_{exact}','FontSize',32,'FontWeight','bold','FontName','Times')
set(get(gca,'Title'),'String','Cyclic Current','FontSize',32,'FontWeight','bold','FontName','Times')
set(gcf,'Position',[1 1 round(1000) round(1000)])

figure
plot(xi,eps_exact(:,Nt*0.25),'--','LineWidth',3); hold on
plot(xi,eps_exact(:,Nt*0.5),'LineWidth',3);
plot(xi,eps_exact(:,Nt*0.75),'LineWidth',3);
plot(xi,eps_exact(:,Nt*1),'-o','LineWidth',3);
legend('\tau=0.25\tau_{total}', '\tau=0.5\tau_{total}', '\tau=0.75\tau_{total}','\tau=\tau_{total}')
axis square
set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'String','\xi','FontSize',42,'FontWeight','bold','FontName','Times')
set(get(gca,'YLabel'),'String','\epsilon_{exact}','FontSize',42,'FontWeight','bold','FontName','Times')
set(get(gca,'Title'),'String','Cyclic Current','FontSize',32,'FontWeight','bold','FontName','Times')
set(gcf,'Position',[1 1 round(1000) round(1000)])



%%========================================================
%% Solution of Error evolution (Paraboloc Equation Finite Difference)

% % compute local residual
% %dI_dtau = omega*t_conv*Im*cos(omega*t_conv.*tau);
% dI_dtau = zeros(Nt,1);
% 
% for i = 1 : size(tau,2) % time
%     for j = 1 : size(xi,2) % space
%         rho(j,i) = dI_dtau(i)*( 0.5*xi(j)^2 - 1/6 - gamma*xi(j)/(1+gamma) + gamma/(2+2*gamma) );
%     end
% end
% 
% figure
% surf(tau, xi,rho,'EdgeColor','none')
% set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
% set(get(gca,'XLabel'),'String','\tau','FontSize',42,'FontWeight','bold','FontName','Times')
% set(get(gca,'YLabel'),'String','\xi','FontSize',42,'FontWeight','bold','FontName','Times')
% set(get(gca,'ZLabel'),'String','\rho','FontSize',42,'FontWeight','bold','FontName','Times')
% set(get(gca,'Title'),'String','Constant Current','FontSize',32,'FontWeight','bold','FontName','Times')
% set(gcf,'Position',[1 1 round(1000) round(1000)])
% 
% figure
% plot(xi,rho(:,Nt*0.25),'--','LineWidth',3); hold on
% plot(xi,rho(:,Nt*0.5),'LineWidth',3);
% plot(xi,rho(:,Nt*0.75),'LineWidth',3);
% plot(xi,rho(:,Nt*1),'-o','LineWidth',3);
% legend('\tau=0.25\tau_{total}', '\tau=0.5\tau_{total}', '\tau=0.75\tau_{total}','\tau=\tau_{total}')
% axis square
% set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
% set(get(gca,'XLabel'),'String','\xi','FontSize',42,'FontWeight','bold','FontName','Times')
% set(get(gca,'YLabel'),'String','{\rho}','FontSize',42,'FontWeight','bold','FontName','Times')
% set(get(gca,'Title'),'String','Constant Current','FontSize',32,'FontWeight','bold','FontName','Times')
% set(gcf,'Position',[1 1 round(1000) round(1000)])
% 
% % % compute IC: epsilon(tau=0) = etaHF(tau=0) - etaLF(tau=0): vector 1*Nx
% etaHF(1:Nx,1) = 0;
% ICond = etaHF(:,1) - etaLF(:,1);
% 
% % % compute BCs: 
% % %-- deta/dxi = alpha  at xi=0
% alpha_error = zeros(1,Nt);
%  
% % %-- deta/dxi = alpha  at xi=1
% beta_error = zeros(1,Nt);
% 
% 
% % compute evolution of error
% epsilon = error_evolution(tau(1),tau(end),xi(1),xi(end), Ntau,Nxi, rho, ICond, alpha_error, beta_error);
% 
% figure
% surf(tau, xi,epsilon,'EdgeColor','none')
% % axis square
% set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
% set(get(gca,'XLabel'),'String','\tau','FontSize',42,'FontWeight','bold','FontName','Times')
% set(get(gca,'YLabel'),'String','\xi','FontSize',42,'FontWeight','bold','FontName','Times')
% set(get(gca,'ZLabel'),'String','\epsilon','FontSize',42,'FontWeight','bold','FontName','Times')
% set(get(gca,'Title'),'String','Constant Current','FontSize',32,'FontWeight','bold','FontName','Times')
% set(gcf,'Position',[1 1 round(1000) round(1000)])
% 
% figure
% plot(xi,epsilon(:,Nt*0.25),'--','LineWidth',3); hold on
% plot(xi,epsilon(:,Nt*0.5),'LineWidth',3);
% plot(xi,epsilon(:,Nt*0.75),'LineWidth',3);
% plot(xi,epsilon(:,Nt*1),'-o','LineWidth',3);
% legend('\tau=0.25\tau_{total}', '\tau=0.5\tau_{total}', '\tau=0.75\tau_{total}','\tau=\tau_{total}')
% axis square
% set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
% set(get(gca,'XLabel'),'String','\xi','FontSize',42,'FontWeight','bold','FontName','Times')
% set(get(gca,'YLabel'),'String','{\epsilon}','FontSize',42,'FontWeight','bold','FontName','Times')
% set(get(gca,'Title'),'String','Constant Current','FontSize',32,'FontWeight','bold','FontName','Times')
% set(gcf,'Position',[1 1 round(1000) round(1000)])


%%========================================================
%% Computing Vcell
etaHF_xi0 = etaHF(1,:);
etaHF_xi1 = etaHF(end,:);
etaLF_xi0 = etaLF(1,:);
etaLF_xi1 = etaLF(end,:);

V_elecHF = ((1+2*gamma)/(1+gamma))*etaHF_xi1 - (gamma/(1+gamma))*etaHF_xi0 - (gamma/(1+gamma)^2)*I;
V_elecLF = ((1+2*gamma)/(1+gamma))*etaLF_xi1 - (gamma/(1+gamma))*etaLF_xi0 - (gamma/(1+gamma)^2)*I;


B = (Ls/kappa_s)*(kappa*sigma/(L*(kappa+sigma)));

V_cellHF = 1 - 0.5*B*I - V_elecHF;
V_cellLF = 1 - 0.5*B*I - V_elecLF;


figure
plot(tau,V_cellHF,'--','LineWidth',3); hold on
plot(tau,V_cellLF,'LineWidth',3); hold on
legend('HF model', 'LF model')
axis square
set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'String','\tau','FontSize',42,'FontWeight','bold','FontName','Times')
set(get(gca,'YLabel'),'String','V^*_{cell}','FontSize',42,'FontWeight','bold','FontName','Times')
set(get(gca,'Title'),'String','Cyclic Current','FontSize',32,'FontWeight','bold','FontName','Times')
set(gcf,'Position',[1 1 round(1000) round(1000)])

