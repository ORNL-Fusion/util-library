clearvars;


suffix = '96'; element='c';
% suffix = '96'; element='h';
fname = ['C:\Work\ADAS\adf11_all\plt',suffix,'\','plt',suffix,'_',element,'.dat'];  % Radiated power (W cm^3)
plt = read_adas_adf11_file(fname);
fname = ['C:\Work\ADAS\adf11_all\prb',suffix,'\','prb',suffix,'_',element,'.dat'];  % Radiated power (W cm^3)
prb = read_adas_adf11_file(fname);



% Te_test = linspace(0.1,5,100);
Te_test = logspace(-1,4,100);  % eV
ne_test = 1e13;   % cm^-3


nCharge = plt.charge;
for iCharge = 1:nCharge
    for i=1:length(ne_test)
        plt_eval(:,iCharge) = interp_adas_rate_coefficient(Te_test,ne_test,plt.te,plt.ne,plt.coeff(:,:,iCharge));
        prb_eval(:,iCharge) = interp_adas_rate_coefficient(Te_test,ne_test,prb.te,prb.ne,prb.coeff(:,:,iCharge));
    end
end


figure; hold on; box on; grid on;
plot(Te_test,plt_eval/1e6)
% plot(Te_test,sum(plt_eval,2)/1e6,'k')
plot(Te_test,prb_eval/1e6)
% plot(Te_test,(prb_eval+prb_eval)/1e6)
% plot(Te_test,sum(prb_eval,2)/1e6,'k')
set(gca,'xscale','log','yscale','log')
 xlabel('T_e (eV)','fontsize',12)
ylabel('<\sigmav> (Wm^3)','fontsize',12)
legend('0','1','2','3','4','5','tot')
axis([0.1,1e4,1e-35,1e-30])

% asdfasd
% s = styflipper(length(ne_test));
% for i=1:length(ne_test)
%     plot(Te_test,sv_iz(:,i),'r','linewidth',2,'linestyle',char(s{i}))
%     plot(Te_test,sv_rc(:,i),'b','linewidth',2,'linestyle',char(s{i}))
% end
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlabel('T_e (eV)','fontsize',12)
% ylabel('<\sigmav> (cm^3/s)','fontsize',12)
% set(gca,'fontsize',12)
% % axis([1,2e4,10^-12,10^-7])
% % axis([0,5,10^-16,10^-8])
% axis([1,2e4,10^-16,10^-7])

