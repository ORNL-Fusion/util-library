clearvars;

GET_COIL_CURRENTS = 1;  % required to get config
GET_ALL_TS = 1;
GET_ITOR = 1;
GET_XICS = 1;

USE_REAL_TS = 1;
PLOT_IT = 1;
SAVE_DATA = 1;
save_dir = '\\sv-it-fs-1\roaming$\jerl\Documents\MATLAB\OP12a_Data';


mydate = '29-Nov-2017';

myshots = [36:44];

% -------------- Get shots and indices -----------------------
fname_shotlist = ['shots','.',mydate,'.mat'];
if ~exist(fname_shotlist,'file')
    fprintf('Creating shotlist matfile: %s\n',fname_shotlist)
    shots = W7X_exp_inq([mydate,' 00:00:00'],[mydate,' 23:59:59'],'quiet');
    save(fname_shotlist,'shots');
else
    fprintf('Reading shotlist matfile: %s\n',fname_shotlist)
    tmp = load(fname_shotlist);
    shots = tmp.shots;
end

for myshot = myshots
    
    clear times currents config ecrh wdia itor prad interf thom ece xics
    
    save_name = fullfile(save_dir,[mydate,'_Shot_',num2str(myshot),'.mat']);
    
    
    % Confirm indexing is sane
    [~,~,c] = fileparts(shots.id{myshot}); c = str2double(c(2:end));
    if ~(c == myshot)
        error('Indexing is insane')
    else
        fprintf('Using shot %d from date %s\n',myshot,mydate)
    end
    
    times.time1 = shots.from(myshot);
    times.time2 = shots.upto(myshot);
    times.t1 = shots.t1(myshot);
    times.t6 = shots.t6(myshot);
    % -------------- \Get shots and indices -----------------------
    
    % ---------------------- Coils and config ----------------------
    if GET_COIL_CURRENTS
        npsc1 = W7X_CoDaC_extract('NPSC1',times.time1,times.time2);
        npsc2 = W7X_CoDaC_extract('NPSC2',times.time1,times.time2);
        npsc3 = W7X_CoDaC_extract('NPSC3',times.time1,times.time2);
        npsc4 = W7X_CoDaC_extract('NPSC4',times.time1,times.time2);
        npsc5 = W7X_CoDaC_extract('NPSC5',times.time1,times.time2);
        psca = W7X_CoDaC_extract('PSCA',times.time1,times.time2);
        pscb = W7X_CoDaC_extract('PSCB',times.time1,times.time2);
        acms11 = W7X_CoDaC_extract('ACM11_ACT',times.time1,times.time2);
        acms19 = W7X_CoDaC_extract('ACM19_ACT',times.time1,times.time2);
        %     acms21 = W7X_CoDaC_extract('ACM21_ACT',times.time1,times.time2);
        %     acms29 = W7X_CoDaC_extract('ACM29_ACT',times.time1,times.time2);
        %     acms31 = W7X_CoDaC_extract('ACM31_ACT',times.time1,times.time2);
        %     acms39 = W7X_CoDaC_extract('ACM39_ACT',times.time1,times.time2);
        %     acms41 = W7X_CoDaC_extract('ACM41_ACT',times.time1,times.time2);
        %     acms49 = W7X_CoDaC_extract('ACM49_ACT',times.time1,times.time2);
        %     acms51 = W7X_CoDaC_extract('ACM51_ACT',times.time1,times.time2);
        %     acms59 = W7X_CoDaC_extract('ACM59_ACT',times.time1,times.time2);
        
        currents.npc1 = round(mean(npsc1.data));
        currents.npc2 = round(mean(npsc2.data));
        currents.npc3 = round(mean(npsc3.data));
        currents.npc4 = round(mean(npsc4.data));
        currents.npc5 = round(mean(npsc5.data));
        currents.pca = round(mean(psca.data));
        currents.pcb = round(mean(pscb.data));
        currents.acm11 = round(mean(acms11.data));
        currents.acm19 = round(mean(acms19.data));
        %     currents.acm21 = round(mean(acms21.data));
        %     currents.acm29 = round(mean(acms29.data));
        %     currents.acm31 = round(mean(acms31.data));
        %     currents.acm39 = round(mean(acms39.data));
        %     currents.acm41 = round(mean(acms41.data));
        %     currents.acm49 = round(mean(acms49.data));
        %     currents.acm51 = round(mean(acms51.data));
        %     currents.acm59 = round(mean(acms59.data));
        
        
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
        currents.taper = [currents.npc1,currents.npc2,currents.npc3,currents.npc4,currents.npc5,currents.pca,currents.pcb,currents.acm11,-currents.acm19];
        
        config=get_w7x_fixedboundary(currents.npc1,currents.npc2,currents.npc3,currents.npc4,currents.npc5,currents.pca,currents.pcb);
        fprintf('\nConfig: %s\n',config.x3_letter_code);
        fprintf('B_0 = %.3f [T]\n',str2double(config.Bax_phi_0__T))
        fprintf('iota_0 = %.3f\n',str2double(config.centralIota))
        fprintf('mirror = %.1f %%\n',str2double(config.mirrorRatio)*100)
        
        if GET_XICS
            FIND_AXIS = 1;
        else
            FIND_AXIS = 0;
        end
        if FIND_AXIS
            coils_file = 'C:\Work\coils.w7x';
            coil = load_vmec_coils_file(coils_file);
            coil = set_w7x_current(coil,currents.taper); % taper = [I1,I2,I3,I4,I5,IA,IB,IS1,IS2];
            bfield.type = 'just_coils';
            bfield.coil = coil.coil;
            bfield.current = coil.current;
            bfield.nsym = coil.num_periods;
            
            phistart = 0*pi/180; Rax_guess = 5.92; Zax_guess = 0;
            [config.Rax,config.Zax] = find_axis(bfield,phistart,Rax_guess,Zax_guess);
            [config.Bax] = bfield_general_rzphi(config.Rax,config.Zax,phistart,bfield);
            config.Bax.B0 = sqrt(config.Bax.br.^2 + config.Bax.bz.^2 + config.Bax.bphi.^2);
            fprintf('B0 = %f\n',config.Bax.B0)
            
            % Make axis fieldline
            nowarn = 0;
            dphi = 0.5*pi/180;
            rnst = 2*pi/bfield.nsym/dphi;
            nsteps = round(rnst);
            config.fl_ax = follow_fieldlines_rzphi_dphi(bfield,config.Rax,config.Zax,phistart,dphi,nsteps,nowarn);
        end
        
    end
    % ---------------------- \Coils and config ----------------------
    
    
    % ------------------ ECRH and plasma time indices --------------
    temp_path='http://archive-webapi.ipp-hgw.mpg.de/ArchiveDB/raw/W7X/CBG_ECRH/TotalPower_DATASTREAM/VX/0/';
    for i=3:-1:1
        temp_path2=strrep(temp_path,'/VX/',['/V' num2str(i,'%1.1i') '/']);
        ecrh = W7X_CoDaC_extract(temp_path2,times.time1,times.time2);
        if ~isempty(ecrh), break; end
    end
    ecrh.totalecrh = mean(ecrh.data(ecrh.data>0.1D3))./1e3; %MW
    dt = diff(double(ecrh.time-times.t1)./1D9);
    ecrh.Winj = sum(ecrh.data(2:end).*dt)./1e3; %MJ
    decrh = diff(ecrh.data);
    temp2 = ecrh.time(decrh >  200);  ecrh.ton  = temp2(1);
    temp2 = ecrh.time(decrh < -200);  ecrh.toff = temp2(end);
    fprintf('Mean Pecrh = %.2f MW, Winj = %.2f MJ \n',ecrh.totalecrh,ecrh.Winj)
    % ------------------ \ECRH and plasma time indices --------------
    
    if 1
        % Wdia
        wdia= W7X_CoDaC_extract('WDIA',times.time1,times.time2); %kJ
        wdia.ecrh_ion =  find(wdia.time>ecrh.ton ,1,'first');
        wdia.ecrh_ioff = find(wdia.time>ecrh.toff,1,'first');
        fprintf('Mean Wdia = %.2f [kJ]\n',mean(wdia.data(wdia.ecrh_ion:wdia.ecrh_ioff)))
        
        if GET_ITOR
            % toroidal current
            itor= W7X_CoDaC_extract('ITOR',times.time1,times.time2); %A
            itor.ecrh_ion =  find(itor.time>ecrh.ton ,1,'first');
            itor.ecrh_ioff = find(itor.time>ecrh.toff,1,'first');
            fprintf('Mean Itor = %.2f [kA] \n',1e-3*mean(itor.data(itor.ecrh_ion:itor.ecrh_ioff)))
        end
        
        % PRAD
        prad=W7X_CoDaC_extract('PRAD',times.time1,times.time2);
        prad.ecrh_ion =  find(prad.time>ecrh.ton ,1,'first');
        prad.ecrh_ioff = find(prad.time>ecrh.toff,1,'first');
        fprintf('Mean Prad = %.2f [MW] \n',1e-6*mean(prad.data(prad.ecrh_ion:prad.ecrh_ioff)))
    end
    
    % %Gas
    % [gasdata,gastime] = get_w7x_gasinlet(times.time1,times.time2,'target');
    
    % Density and temperature
    
    % -- Interferometer
    interf.Ltot = 1.33;
    [interf.data, interf.time] = get_w7x_interferrometer(times.time1,times.time2);
    interf.ecrh_ion =  find(interf.time>ecrh.ton ,1,'first');
    interf.ecrh_ioff = find(interf.time>ecrh.toff,1,'first');
    fprintf('Mean ne (interf) = %.2f [10^19 m^-3] \n',1e-19./interf.Ltot*mean(interf.data(interf.ecrh_ion:interf.ecrh_ioff)))
    
    % Thomson central ne
    temp_path='http://archive-webapi.ipp-hgw.mpg.de/Test/raw/W7X/QTB_Central/volume_2_DATASTREAM/VX/1/ne_map/_signal.html';
    thom = [];
    for i=3:-1:1
        temp_path2=strrep(temp_path,'/VX/',['/V' num2str(i,'%1.1i') '/']);
        thom.central_ne = W7X_CoDaC_extract(temp_path2,times.time1,times.time2);
        if ~isempty(thom.central_ne), break; end
    end
    if ~isempty(thom.central_ne)
        thom.central_ne.ecrh_ion =  find(thom.central_ne.time>ecrh.ton ,1,'first');
        thom.central_ne.ecrh_ioff = find(thom.central_ne.time>ecrh.toff,1,'first');
        
        % Thomson central Te
        temp_path='http://archive-webapi.ipp-hgw.mpg.de/Test/raw/W7X/QTB_Central/volume_2_DATASTREAM/VX/0/Te_map/_signal.html';
        for i=3:-1:1
            temp_path2=strrep(temp_path,'/VX/',['/V' num2str(i,'%1.1i') '/']);
            thom.central_Te = W7X_CoDaC_extract(temp_path2,times.time1,times.time2);
            if ~isempty(thom.central_Te), break; end
        end
        
        if GET_ALL_TS
            fprintf('Getting all TS data for this shot, slowly\n')
            tic1=tic; [thom.time, thom.te, thom.sigma_te, thom.Ne, thom.sigma_Ne, thom.x,thom.y,thom.z] = get_w7x_thomson(times.time1,times.time2); toc(tic1)
            thom.r = sqrt(thom.x.*thom.x+thom.y.*thom.y);
            thom.PHI_THOM_USE = mod(atan2(thom.y,thom.x),2*pi/5);
            thom.zax = interp1(config.fl_ax.phi,config.fl_ax.z,thom.PHI_THOM_USE);
            thom.rax = interp1(config.fl_ax.phi,config.fl_ax.r,thom.PHI_THOM_USE);
            [~,thom.IND_THOM_USE] = min(sqrt((thom.zax-thom.z).^2 + (thom.rax-thom.r).^2));
            
        end
        
        
        
    end
    
    
    % %ECE
    ece = [];
    %[ece.time, ece.te, ece.sigma_te, ece.freq, ece.bw] = get_w7x_ece(times.time1,times.time2); % eV
    
    
    if GET_XICS
        % XICS
        temp = strsplit(shots.id{myshot},'.');
        num = sscanf([temp{1}(3:end) temp{2}],'%i');
        xics=get_w7x_xics(num);
        phi_use_xics = mod(mean([xics.p1;xics.p2]),2*pi/5);
        xics.rax = interp1(config.fl_ax.phi,config.fl_ax.r,phi_use_xics);
        xics.zax = interp1(config.fl_ax.phi,config.fl_ax.z,phi_use_xics);
        for i = 1:length(xics.r1)
            mm = (xics.z2(i)-xics.z1(i))/(xics.r2(i)-xics.r1(i));
            bint = xics.z1(i)-mm*xics.r1(i);
            aa = mm;
            bb = -1;
            cc = bint;
            dxics(i) = abs(aa*xics.rax + bb*xics.zax + cc)/sqrt(aa^2+bb^2);
        end
        [~,xics.ch_use] = min(dxics);
    end
    
    if SAVE_DATA
        save(save_name,'myshot','times','currents','config','ecrh',...
            'wdia','itor','prad','interf','thom','ece','xics');
        fprintf('Saved %s\n',save_name)
    end
    
    if PLOT_IT
        nsub = 3; fs = 14;
        %         if myshot == myshots(1)
        figure; hold on; box on;
        %
        %         end
        subplot(nsub,1,1); hold on; box on; set(gca,'fontsize',fs)
        title(sprintf('Shot = %i',myshot))
        set(gca,'xlim',[-0.1,double(ecrh.toff-ecrh.ton)*1e-9+1])
        plot(double(ecrh.time - ecrh.ton)*1e-9,ecrh.data/1e3)
        plot(double(prad.time  - ecrh.ton)*1e-9,prad.data/1e6)
        yy = get(gca,'ylim'); set(gca,'ylim',[0,yy(2)]);
        ylabel('Power (MW)','fontsize',fs)
        set(gca,'xlim',[-0.1,double(ecrh.toff-ecrh.ton)*1e-9+1])
        
        subplot(nsub,1,2); hold on; box on; set(gca,'fontsize',fs)
        set(gca,'xlim',[-0.1,double(ecrh.toff-ecrh.ton)*1e-9+1])
        plot(double(interf.time-ecrh.ton)*1e-9,interf.data./interf.Ltot*1e-19)
        if ~isempty(thom.central_ne)
            plot(double(thom.central_ne.time-ecrh.ton)*1e-9,thom.central_ne.data,'.-')
            if USE_REAL_TS
                plot(double(thom.time-ecrh.ton)*1e-9,thom.Ne(:,thom.IND_THOM_USE),'x-')
            end
        end
        ylabel('n_e (10^{13} m^{-3})','fontsize',fs)
        set(gca,'xlim',[-0.1,double(ecrh.toff-ecrh.ton)*1e-9+1])
        
        subplot(nsub,1,3); hold on; box on; set(gca,'fontsize',fs)
        set(gca,'xlim',[-0.1,double(ecrh.toff-ecrh.ton)*1e-9+1])
        if ~isempty(thom.central_ne)
            tuse = double(thom.central_Te.time-ecrh.ton)*1e-9;
            
            tuse_inds = tuse <(double(ecrh.toff-ecrh.ton)*1e-9+0.5);
            plot(tuse(tuse_inds),thom.central_Te.data(tuse_inds),'.-')
            
            if USE_REAL_TS
                plot(double(thom.time-ecrh.ton)*1e-9,thom.te(:,thom.IND_THOM_USE),'x-')
            end
        end
        %plot(double(ece.time - ecrh.ton)*1e-9,ece.te)
        % ylabel('T (keV)','fontsize',fs)
        if GET_XICS
            plot(xics.ti_time,xics.ti(:,xics.ch_use))
        end
        % plot(xics.te_time,xics.te(:,xics_ch_use))
        yy = get(gca,'ylim'); set(gca,'ylim',[0,min(yy(2),10)]);
        
        
        % subplot(nsub,1,4); hold on; box on; set(gca,'fontsize',fs)
        % set(gca,'xlim',[-0.1,double(ecrh.toff-ecrh.ton)*1e-9+1])
        plot(double(wdia.time-ecrh.ton)*1e-9,wdia.data./100,'-')
        xlabel('t - t_{ECRH}^{on} (s)','fontsize',fs)
        % ylabel('W_{dia} (kJ)','fontsize',fs)
        drawnow;
    end
end
