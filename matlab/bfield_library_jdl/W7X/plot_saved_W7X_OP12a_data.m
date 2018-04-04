clearvars;

save_dir = 'C:\Users\jjl\Dropbox (ORNL)\W7X\OP12a_Data';

mydate = '29-Nov-2017';
  
% myshots = [41,39,32,44,21];  % H 2MW density scan (0,11,22,32,43)
% myshots = [41,42];

% myshots = [40,37,27,43,10]; % He 2MW density scan (0,11,22,32,43)

myshots = [36,30]
JUST_BASICS = 0;
USE_REAL_TS = 1;
GET_XICS = 0;

coils_file = 'C:\Work\Stellarator\W7X EMC3 modeling\Mark coils and input\coils.w7x';
coil = load_vmec_coils_file(coils_file);


for myshot = myshots
    
    clear times currents config ecrh wdia itor prad interf thom ece xics
    
    save_name = fullfile(save_dir,[mydate,'_Shot_',num2str(myshot),'.mat']);
    load(save_name);
    
    fprintf('Coil currents ---------------------\n')
    fprintf('I1 = %.0f\n',currents.npc1)
    fprintf('I2 = %.0f\n',currents.npc2)
    fprintf('I3 = %.0f\n',currents.npc3)
    fprintf('I4 = %.0f\n',currents.npc4)
    fprintf('I5 = %.0f\n',currents.npc5)
    fprintf('IA = %.0f\n',currents.pca)
    fprintf('IB = %.0f\n',currents.pcb)
    fprintf('IS1 = %.0f\n',currents.acm11)
    fprintf('IS2 = %.0f\n',-currents.acm19)
    fprintf('IS COIL SIGNS ONLY TO GUESSED TO VMEC CONVENTION!!!\n')
    %         currents.taper = [currents.npc1,currents.npc2,currents.npc3,currents.npc4,currents.npc5,currents.pca,currents.pcb,currents.acm11,-currents.acm19];
    
    
    fprintf('\nConfig: %s\n',config.x3_letter_code);
    fprintf('B_0 = %.3f [T]\n',str2double(config.Bax_phi_0__T))
    fprintf('iota_0 = %.3f\n',str2double(config.centralIota))
    fprintf('mirror = %.1f %%\n',str2double(config.mirrorRatio)*100)
    
    
    coil = set_w7x_current(coil,currents.taper); % taper = [I1,I2,I3,I4,I5,IA,IB,IS1,IS2];
    bfield.type = 'just_coils';
    bfield.coil = coil.coil;
    bfield.current = coil.current;
    bfield.nsym = coil.num_periods;
    
    
    fprintf('Mean Pecrh = %.2f MW, Winj = %.2f MJ \n',ecrh.totalecrh,ecrh.Winj)
    
    fprintf('Mean Wdia = %.2f [kJ]\n',mean(wdia.data(wdia.ecrh_ion:wdia.ecrh_ioff)))
    
    fprintf('Mean Itor = %.2f [kA] \n',1e-3*mean(itor.data(itor.ecrh_ion:itor.ecrh_ioff)))
    
    fprintf('Mean Prad = %.2f [MW] \n',1e-6*mean(prad.data(prad.ecrh_ion:prad.ecrh_ioff)))
    
    fprintf('Mean ne (interf) = %.2f [10^19 m^-3] \n',1e-19./interf.Ltot*mean(interf.data(interf.ecrh_ion:interf.ecrh_ioff)))
    
    nsub = 3; fs = 14;
%             if myshot == myshots(1)
    figure; hold on; box on;
    %
%             end
    subplot(nsub,1,1); hold on; box on; set(gca,'fontsize',fs)
    title(sprintf('Shot = %i',myshot))
    set(gca,'xlim',[-0.1,double(ecrh.toff-ecrh.ton)*1e-9+1])
    plot(double(ecrh.time - ecrh.ton)*1e-9,ecrh.data/1e3)
    if ~JUST_BASICS
    plot(double(prad.time  - ecrh.ton)*1e-9,prad.data/1e6)
    legend('ECRH','PRAD')
    end
    yy = get(gca,'ylim'); set(gca,'ylim',[0,yy(2)]);
    ylabel('Power (MW)','fontsize',fs)
    set(gca,'xlim',[-0.1,double(ecrh.toff-ecrh.ton)*1e-9+1])
    
    
    subplot(nsub,1,2); hold on; box on; set(gca,'fontsize',fs)
    set(gca,'xlim',[-0.1,double(ecrh.toff-ecrh.ton)*1e-9+1])
    plot(double(interf.time-ecrh.ton)*1e-9,interf.data./interf.Ltot*1e-19)
    if ~JUST_BASICS
    if ~isempty(thom.central_ne)
        plot(double(thom.central_ne.time-ecrh.ton)*1e-9,thom.central_ne.data,'.-')
        if USE_REAL_TS
            plot(double(thom.time-ecrh.ton)*1e-9,thom.Ne(:,thom.IND_THOM_USE),'x-')
        end
    end
    legend('interf','TS uncal?','TS cal?')
    end
    ylabel('n_e (10^{13} m^{-3})','fontsize',fs)
    set(gca,'xlim',[-0.1,double(ecrh.toff-ecrh.ton)*1e-9+1])
    
    
    subplot(nsub,1,3); hold on; box on; set(gca,'fontsize',fs)
    set(gca,'xlim',[-0.1,double(ecrh.toff-ecrh.ton)*1e-9+1])
    plot(double(wdia.time-ecrh.ton)*1e-9,wdia.data./100,'-')
    if ~JUST_BASICS
    if ~isempty(thom.central_ne)
        tuse = double(thom.central_Te.time-ecrh.ton)*1e-9;
        
        tuse_inds = tuse <(double(ecrh.toff-ecrh.ton)*1e-9+0.5);
        plot(tuse(tuse_inds),thom.central_Te.data(tuse_inds),'.-')
        
        if USE_REAL_TS
            plot(double(thom.time-ecrh.ton)*1e-9,thom.te(:,thom.IND_THOM_USE),'x-')
        end
        legend('Wdia','TS','TS cal?')
    end
    %plot(double(ece.time - ecrh.ton)*1e-9,ece.te)
    % ylabel('T (keV)','fontsize',fs)
    if GET_XICS
        plot(xics.ti_time,xics.ti(:,xics.ch_use))
    end
    % plot(xics.te_time,xics.te(:,xics_ch_use))
    end
    yy = get(gca,'ylim'); set(gca,'ylim',[0,min(yy(2),10)]);
    
    
    % subplot(nsub,1,4); hold on; box on; set(gca,'fontsize',fs)
    % set(gca,'xlim',[-0.1,double(ecrh.toff-ecrh.ton)*1e-9+1])
    
    xlabel('t - t_{ECRH}^{on} (s)','fontsize',fs)
    % ylabel('W_{dia} (kJ)','fontsize',fs)
    
    drawnow;
end
% legend