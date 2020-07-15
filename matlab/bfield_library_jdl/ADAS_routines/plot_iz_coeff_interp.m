clearvars;

% suffix = '12'; element='h';
suffix = '96'; element='he';
fname = ['C:\Work\ADAS\adf11_all\scd',suffix,'\','scd',suffix,'_',element,'.dat'];  % Effective ionization coefficients (cm^-3/s)
scd = read_adas_adf11_file(fname); 

scd1 = squeeze(scd.coeff(:,:,1));  % Choose charge 

Te_test = linspace(.1,20,100);
ne_test = logspace(log10(1e11),log10(1e13),3);  % cm3

for i=1:length(ne_test)
    sv_iz_adas(:,i) = interp_adas_rate_coefficient(Te_test,ne_test(i),scd.te,scd.ne,scd1);
end


s = styflipper(length(ne_test));

figure; hold on; box on
for i=1:length(ne_test)
    plot(Te_test,sv_iz_adas(:,i)./1e6,'b','linewidth',2,'linestyle',char(s{i}))
end
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('T_e (eV)','fontsize',12)
ylabel('<\sigmav> (m^3/s)','fontsize',12)
set(gca,'fontsize',12)
axis([.1,10,10^-20,10^-14])
legend('n_e = 1e17','n_e = 1e18','n_e = 1e19')