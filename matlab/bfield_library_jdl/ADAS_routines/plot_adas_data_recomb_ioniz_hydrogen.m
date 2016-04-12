clearvars;

% suffix = '96'; element='h';
suffix = '50'; element='w';

fprintf('Reading scd\n');
fname = ['C:\Work\ADAS\adf11_all\scd',suffix,'\','scd',suffix,'_',element,'.dat'];  % Effective ionization coefficients (cm^-3/s)
[te_scd,ne_scd,scd] = read_adas_adf11_file(fname);
fprintf('Reading acd\n');
fname = ['C:\Work\ADAS\adf11_all\acd',suffix,'\','acd',suffix,'_',element,'.dat'];  % Effective recombination coefficients (cm^-3/s)
[te_acd,ne_acd,acd] = read_adas_adf11_file(fname);
fprintf('Reading plt\n');
fname = ['C:\Work\ADAS\adf11_all\plt',suffix,'\','plt',suffix,'_',element,'.dat'];  % Radiated power (W cm^3)
[te_plt,ne_plt,plt] = read_adas_adf11_file(fname);

ne_plot_targ = [1e13,1e14,1e15];  % cm^-3

for i = 1:length(ne_plot_targ)
    [~,in_scd]=min(abs(ne_scd - ne_plot_targ(i)));
    [~,in_acd]=min(abs(ne_acd - ne_plot_targ(i)));
    [~,in_plt]=min(abs(ne_plt - ne_plot_targ(i)));
    scd_arr(:,i) = squeeze(scd(:,in_scd,:));
    acd_arr(:,i) = squeeze(acd(:,in_acd,:));
    fprintf('Plotting scd for ne = %e\n',ne_scd(in_scd))
    fprintf('Plotting acd for ne = %e\n',ne_acd(in_acd))
    fprintf('Plotting plt for ne = %e\n',ne_plt(in_plt))
end


s = styflipper(length(ne_plot_targ));
figure; hold on; box on;
for i = 1:length(ne_plot_targ)
    plot(te_scd,scd_arr(:,i),'r','linewidth',2,'linestyle',char(s{i}))
    plot(te_acd,acd_arr(:,i),'b','linewidth',2,'linestyle',char(s{i}))
end
% plot(te_acd,squeeze(acd(1,:,:)),'--')
% set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('T_e (eV)','fontsize',12)
ylabel('<\sigmav> (cm^3/s)','fontsize',12)
set(gca,'fontsize',12)
axis([0,5,10^-16,10^-8])
legend('iz','di-el rc')
% set(gca,'linewidth',2)



% figure; hold on; box on;
% plot(te_plt,squeeze(plt(:,in_plt,:)),'-')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlabel('T_e (eV)','fontsize',12)
% ylabel('<\sigmav> (cm^3/s)','fontsize',12)
% set(gca,'fontsize',12)
% % axis([1,1e4,10^-12,10^-6])
% 


