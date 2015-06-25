clearvars;

% cc = cc(iz,ne,te)
% fname = ['C:\Work\ADAS\adf11_all\scd96\','scd96_h.dat'];  % Effective ionization coefficients (cm^-3/s)
% [te_scd,ne_scd,scd] = read_adas_adf11_file(fname);
% fname = ['C:\Work\ADAS\adf11_all\acd96\','acd96_h.dat'];  % Effective recombination coefficients (cm^-3/s)
% [te_acd,ne_acd,acd] = read_adas_adf11_file(fname);
% fname = ['C:\Work\ADAS\adf11_all\plt96\','plt96_h.dat'];  % Radiated power (W cm^3)
% [te_plt,ne_plt,plt] = read_adas_adf11_file(fname);



fname = ['C:\Work\ADAS\adf11_all\scd96\','scd96_n.dat'];  % Effective ionization coefficients (cm^-3/s)
[te_scd,ne_scd,scd] = read_adas_adf11_file(fname);
fname = ['C:\Work\ADAS\adf11_all\acd96\','acd96_n.dat'];  % Effective recombination coefficients (cm^-3/s)
[te_acd,ne_acd,acd] = read_adas_adf11_file(fname);
% fname = ['C:\Work\ADAS\adf11_all\plt96\','plt96_n.dat'];  % Radiated power (W cm^3)
fname = ['C:\Work\ADAS\','plt96_c.dat'];  % Radiated power (W cm^3)
[te_plt,ne_plt,plt] = read_adas_adf11_file(fname);

ne_plot_targ = 1e13;  % cm^-3
[~,in_scd]=min(abs(ne_scd - ne_plot_targ));
[~,in_acd]=min(abs(ne_acd - ne_plot_targ));
[~,in_plt]=min(abs(ne_plt - ne_plot_targ));
fprintf('Plotting scd for ne = %e\n',ne_scd(in_scd))
fprintf('Plotting acd for ne = %e\n',ne_acd(in_acd))
fprintf('Plotting plt for ne = %e\n',ne_plt(in_plt))

figure; hold on; box on;
plot(te_scd,squeeze(scd(:,in_scd,:)),'-')
plot(te_acd,squeeze(acd(:,in_acd,:)),'--')
% plot(te_acd,squeeze(acd(1,:,:)),'--')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('T_e (eV)','fontsize',12)
ylabel('<\sigmav> (cm^3/s)','fontsize',12)
set(gca,'fontsize',12)
axis([1,1e4,10^-12,10^-6])


figure; hold on; box on;
plot(te_plt,squeeze(plt(:,in_plt,:)),'-')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('T_e (eV)','fontsize',12)
ylabel('<\sigmav> (cm^3/s)','fontsize',12)
set(gca,'fontsize',12)
% axis([1,1e4,10^-12,10^-6])



