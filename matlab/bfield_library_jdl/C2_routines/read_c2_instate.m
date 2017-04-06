% function read_c2_instate(fname)
clearvars;
fname = 'C:\Work\SCREAM\FASTRAN\FASTRAN_input\';

icheck = exist(fname);
if icheck == 0
    error('File/dir does not exist: %s\n',fname)
end
if icheck == 7
    fname = fullfile(fname,'instate');   
end

instate=read_namelist(fname,'INSTATE');

figure; hold on; box on;
subplot(3,1,1); hold on; box on;
plot(instate.rho,instate.ne,'linewidth',3)
ylabel('n_e [10^{19} m^{-3}]','fontsize',14,'fontweight','bold')
set(gca,'fontsize',14,'fontweight','bold')
subplot(3,1,2); hold on; box on;
plot(instate.rho,instate.te,'linewidth',3)
plot(instate.rho,instate.ti,'linewidth',3)
ylabel('T_i [eV]','fontsize',14,'fontweight','bold')
set(gca,'fontsize',14,'fontweight','bold')
subplot(3,1,3); hold on; box on;
plot(instate.rho,instate.omega/1e3,'linewidth',3)
xlabel('\rho','fontsize',14,'fontweight','bold')
set(gca,'fontsize',14,'fontweight','bold')
ylabel('\omega [krad/s]','fontsize',14,'fontweight','bold')