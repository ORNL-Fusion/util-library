% test_adas_interp
clearvars;

fname = ['C:\Work\ADAS\adf11_all\scd96\','scd96_h.dat'];  % Effective ionization coefficients (cm^-3/s) 
[te_scd,ne_scd,scd] = read_adas_adf11_file(fname);
fname = ['C:\Work\ADAS\adf11_all\acd96\','acd96_h.dat'];  % Effective recombination coefficients (cm^-3/s)
[te_acd,ne_acd,acd] = read_adas_adf11_file(fname);
fname = ['C:\Work\ADAS\adf11_all\ccd96\','ccd96_h.dat'];  % Effective cx coefficients (cm^-3/s)
[te_ccd,ne_ccd,ccd] = read_adas_adf11_file(fname);


% Te_test = linspace(0.1,5,100);
Te_test = logspace(-1,4,100);
ne_test = logspace(13,15,3);

scd1 = squeeze(scd(1,:,:));
acd1 = squeeze(acd(1,:,:));
ccd1 = squeeze(ccd(1,:,:));
for i=1:length(ne_test)
    sv_iz(:,i) = interp_adas_rate_coefficient(Te_test,ne_test(i),te_scd,ne_scd,scd1);
    sv_rc(:,i) = interp_adas_rate_coefficient(Te_test,ne_test(i),te_acd,ne_acd,acd1);
    sv_cx(:,i) = interp_adas_rate_coefficient(Te_test,ne_test(i),te_ccd,ne_ccd,ccd1);
end

s = styflipper(length(ne_test));
figure; hold on; box on
for i=1:length(ne_test)
    plot(Te_test,sv_iz(:,i),'r','linewidth',2,'linestyle',char(s{i}))
    plot(Te_test,sv_rc(:,i),'b','linewidth',2,'linestyle',char(s{i}))
%     plot(Te_test,sv_cx(:,i),'g','linewidth',2,'linestyle',char(s{i}))  % This one is vs Ti
end
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('T_e (eV)','fontsize',12)
ylabel('<\sigmav> (cm^3/s)','fontsize',12)
set(gca,'fontsize',12)
% axis([1,2e4,10^-12,10^-7])
% axis([0,5,10^-16,10^-8])
axis([1,2e4,10^-16,10^-7])

