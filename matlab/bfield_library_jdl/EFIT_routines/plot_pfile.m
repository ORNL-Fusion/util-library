function plot_pfile(p)

figure; hold on; box on; set(gcf,'color','w');set(gca,'fontsize',14,'fontweight','bold'); grid on;
plot(p.profiles.ne_psin,p.ne)
plot(p.te_psin,p.te)



