clearvars;

% suffix = '96'; element='h';
suffix = '96'; element='c';
% suffix = '50'; element='w';

fprintf('Reading scd\n');
fname = ['C:\Work\ADAS\adf11_all\scd',suffix,'\','scd',suffix,'_',element,'.dat'];  % Effective ionization coefficients (cm^-3/s)
scd = read_adas_adf11_file(fname);
fprintf('Reading acd\n');
fname = ['C:\Work\ADAS\adf11_all\acd',suffix,'\','acd',suffix,'_',element,'.dat'];  % Effective recombination coefficients (cm^-3/s)
acd = read_adas_adf11_file(fname);
fprintf('Reading plt\n');
fname = ['C:\Work\ADAS\adf11_all\plt',suffix,'\','plt',suffix,'_',element,'.dat'];  % Radiated power (W cm^3)
plt = read_adas_adf11_file(fname);

ne_plot_targ = [1e13,1e14,1e15];  % cm^-3

for i = 1:length(ne_plot_targ)
    [~,in_scd]=min(abs(scd.ne - ne_plot_targ(i)));
    [~,in_acd]=min(abs(acd.ne - ne_plot_targ(i)));
    [~,in_plt]=min(abs(plt.ne - ne_plot_targ(i)));
    scd_arr(:,i) = squeeze(scd.coeff(in_scd,:)); 
    acd_arr(:,i) = squeeze(acd.coeff(in_acd,:));
    fprintf('Plotting scd for ne = %e\n',scd.ne(in_scd))
    fprintf('Plotting acd for ne = %e\n',acd.ne(in_acd))
    fprintf('Plotting plt for ne = %e\n',plt.ne(in_plt))
end


s = styflipper(length(ne_plot_targ));
figure; hold on; box on;
for i = 1:length(ne_plot_targ)
    plot(scd.te,scd_arr(:,i),'r','linewidth',2,'linestyle',char(s{i}))
    plot(acd.te,acd_arr(:,i),'b','linewidth',2,'linestyle',char(s{i}))
end
% plot(acd.te,squeeze(acd(1,:,:)),'--')
% set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('T_e (eV)','fontsize',12)
ylabel('<\sigmav> (cm^3/s)','fontsize',12)
set(gca,'fontsize',12)
axis([0,5,10^-16,10^-8])
legend('iz','di-el rc')
% set(gca,'linewidth',2)



% figure; hold on; box on;
% plot(plt.te,squeeze(plt.coeff(:,in_plt,:)),'-')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlabel('T_e (eV)','fontsize',12)
% ylabel('<\sigmav> (cm^3/s)','fontsize',12)
% set(gca,'fontsize',12)
% % axis([1,1e4,10^-12,10^-6])
% 


