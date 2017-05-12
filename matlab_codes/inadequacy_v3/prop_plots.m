%% Plots properties

set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
set(get(gca,'XLabel'),'FontSize',32,'FontWeight','bold','FontName','Times')
set(get(gca,'YLabel'),'FontSize',32,'FontWeight','bold','FontName','Times')
set(gcf,'Position',[1 1 round(1000) round(1000)])
set(get(gca,'Title'),'FontSize',32,'FontWeight','bold','FontName','Times')
pbaspect([2 1 1])