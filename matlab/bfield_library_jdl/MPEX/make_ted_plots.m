clearvars;

bfield_method = 1;  % 1 is Biot-Savart, 2 is exact



debug_plots = 0;

current_A = 3300;
% current_A = 6400;
current_B = 0;
% config = 'flat';
% helicon_current = 210;

% current_A = 6607;
% current_B = 6607;
% config = 'focus';
% helicon_current = -80;


skimmer = 1;
geo = get_Proto_geometry(0,0,skimmer);
% configs = [{'flat'},{'standard'},{'focus'}];
% configs = [{'flat'}];
configs = [{'focus'}];

verbose = 0;
num_test = 5;
helicon_currents = linspace(-100,100,num_test);

doLCFS = 1;


% MIRROR RATIOS
for iconfig = 1:length(configs)
    config = configs{iconfig};
    for ih = 1:num_test
        if debug_plots
            geo = get_Proto_geometry(1,1,skimmer);
        end
        fprintf('ih %d of %d\n',ih,num_test)
        helicon_current = helicon_currents(ih);
        if bfield_method == 1
            [coil,current] = build_Proto_coils(helicon_current,current_A,current_B,config,verbose);
            bfield.coil = coil;
            bfield.current = current;
            bfield.type = 'just_coils';
        elseif bfield_method == 2
            [coil,current] = build_Proto_coils_jackson(helicon_current,current_A,current_B,config);
            bfield.coil = coil;
            bfield.current = current;
            bfield.type = 'MPEX';
        else
            error('bad bfield_method')
        end
        ratios = eval_mirror_ratios(bfield,geo);
        Ru1s(iconfig,ih) = ratios.R_u1;
        Ru2s(iconfig,ih) = ratios.R_u2;
        Rd5s(iconfig,ih) = ratios.R_d5;
        Rd6s(iconfig,ih) = ratios.R_d6;
        
        %Skimmer
        Rstart = geo.skimmer.ID/2;
        Zstart = (geo.skimmer.z1+geo.skimmer.z2)/2;
        phistart = 0;
        dz = 0.01;
        L = geo.target.z - Zstart;
        nsteps = round(abs(L/dz));
        dz = L/nsteps;
        f = follow_fieldlines_rzphi_dz(bfield,Rstart,Zstart,phistart,dz,nsteps);
        r_skimmer(iconfig,ih) = f.r(end);
        if debug_plots
            plot(f.z,f.r)
        end
        
        %helicon u
        Rstart = geo.helicon.r;
        Zstart = geo.helicon.z1;
        phistart = 0;
        dz = 0.01;
        L = geo.target.z - Zstart;
        nsteps = round(abs(L/dz));
        dz = L/nsteps;
        f = follow_fieldlines_rzphi_dz(bfield,Rstart,Zstart,phistart,dz,nsteps);
        r_u(iconfig,ih) = f.r(end);        
        f = clip_fl_at_vessel(f,geo);
        if isnan(f.r(end))
            r_u_clip(iconfig,ih) = 1;
        else
            r_u_clip(iconfig,ih) = 0;
        end
        if debug_plots
            plot(f.z,f.r)
        end
        
        
        %helicon m
        Rstart = geo.helicon.r;
        Zstart = (geo.helicon.z2-geo.helicon.z1)/2 + geo.helicon.z1;
        phistart = 0;
        dz = 0.01;
        L = geo.target.z - Zstart;
        nsteps = round(abs(L/dz));
        dz = L/nsteps;
        f = follow_fieldlines_rzphi_dz(bfield,Rstart,Zstart,phistart,dz,nsteps);
        r_m(iconfig,ih) = f.r(end);
        f = clip_fl_at_vessel(f,geo);        
        if isnan(f.r(end))
            r_m_clip(iconfig,ih) = 1;
        else
            r_m_clip(iconfig,ih) = 0;
        end        
        if debug_plots
            plot(f.z,f.r)
        end
        
        
        %helicon d
        Rstart = geo.helicon.r;
        Zstart = geo.helicon.z2;
        phistart = 0;
        dz = 0.01;
        L = geo.target.z - Zstart;
        nsteps = round(abs(L/dz));
        dz = L/nsteps;
        f = follow_fieldlines_rzphi_dz(bfield,Rstart,Zstart,phistart,dz,nsteps);
        r_d(iconfig,ih) = f.r(end);
        f = clip_fl_at_vessel(f,geo);        
        if isnan(f.r(end))
            r_d_clip(iconfig,ih) = 1;
        else
            r_d_clip(iconfig,ih) = 0;
        end             
        if debug_plots
            plot(f.z,f.r)
        end
        
        
        %LCFS
        if doLCFS
            shot_tmp{1} = helicon_current;
            shot_tmp{2} = current_A;
            shot_tmp{3} = current_B;
            shot_tmp{4} = config;
            shot_tmp{5} = skimmer;
            f = find_lcfs(shot_tmp);
            r_lcfs(iconfig,ih) = f.r(1);
            if debug_plots
                plot(f.z,f.r)
            end

            %         r_lcfs_skimmer(iconfig,ih) = f.r(1);
            
            
            %         shot_tmp{5} = 0;
            %         flcfs = find_lcfs(shot_tmp);
            %         r_lcfs_noskimmer(iconfig,ih) = f.r(1);
        end
    end
end

% for iconfig = 1:length(configs)
%     config = configs{iconfig};
%     figure; hold on; box on;
%     plot(helicon_currents,1./Ru1s(iconfig,:),'b','linewidth',2)
%     plot(helicon_currents,1./Ru2s(iconfig,:),'m','linewidth',2)
%     plot(helicon_currents,1./Rd5s(iconfig,:),'c','linewidth',2)
%     plot(helicon_currents,1./Rd6s(iconfig,:),'r','linewidth',2)
%     xlabel('I_{helicon} [Amps]','fontsize',14)
%     h=legend('R_{u1}','R_{u2}','R_{d5}','R_{d6}');
%     ylabel('Mirror ratios','fontsize',14)
%     set(gca,'fontsize',14)
%     set(h,'fontsize',14)
%     title(config,'fontsize',14)
% end


% target radial positions
for iconfig = 1:length(configs)
    config = configs{iconfig};
    figure; hold on; box on;
    if skimmer
        plot(helicon_currents,r_skimmer(iconfig,:),'k-','linewidth',2)
    end        
    if doLCFS
        plot(helicon_currents,r_lcfs(iconfig,:),'m-','linewidth',2)
    end        
    inc = find(r_d_clip(iconfig,:) ~= 1);
    plot(helicon_currents(inc),r_d(iconfig,inc),'r-','linewidth',2)
    inc = find(r_m_clip(iconfig,:) ~= 1);
    plot(helicon_currents(inc),r_m(iconfig,inc),'g-','linewidth',2)
    inc = find(r_u_clip(iconfig,:) ~= 1);    
    if ~isempty(inc)
        plot(helicon_currents(inc),r_u(iconfig,inc),'b-','linewidth',2)
    else
        plot([0,0],[0,0],'b-','linewidth',2)
    end
    ic = find(r_d_clip(iconfig,:) == 1);
    if ic(1) > 1
        ic = [ic(1)-1,ic];
    end    
    plot(helicon_currents(ic),r_d(iconfig,ic),'r--','linewidth',2)
    ic = find(r_m_clip(iconfig,:) == 1);             
    if ic(1) > 1
        ic = [ic(1)-1,ic];
    end    
    plot(helicon_currents(ic),r_m(iconfig,ic),'g--','linewidth',2)
    ic = find(r_u_clip(iconfig,:) == 1);
    if ic(1) > 1
        ic = [ic(1)-1,ic];
    end
    plot(helicon_currents(ic),r_u(iconfig,ic),'b--','linewidth',2)
    
%     plot(helicon_currents,r_m(iconfig,:),'g','linewidth',2)
%     plot(helicon_currents,r_u(iconfig,:),'b','linewidth',2)

    xlabel('I_{helicon} [Amps]','fontsize',14)
    leg_str = [{'r_d'},{'r_m'},{'r_u'}];
    if doLCFS
        leg_str = [{'r_{lcfs}'},leg_str];
    end
    if skimmer
        leg_str = [{'r_{skimmer}'},leg_str];
    end        
    h=legend(leg_str);
    ylabel('r_{target} [m]','fontsize',14)
    set(gca,'fontsize',14)
    set(h,'fontsize',14)
    title(['Target 7.5, ',config],'fontsize',14)
    axis([min(helicon_currents),max(helicon_currents),0,0.05])
end

colors = ['b','r','k'];

% figure; hold on; box on;
% for iconfig=1:length(configs)
%     plot(helicon_currents,1./Ru1s(iconfig,:),colors(iconfig),'linewidth',2)
% end
% xlabel('I_{helicon} [Amps]','fontsize',14)
% h=legend(configs);
% set(gca,'fontsize',14)
% set(h,'fontsize',14)
% ylabel('R_{u1}','fontsize',14)
% 
% figure; hold on; box on;
% for iconfig=1:length(configs)
%     plot(helicon_currents,1./Rd6s(iconfig,:),colors(iconfig),'linewidth',2)
% end
% xlabel('I_{helicon} [Amps]','fontsize',14)
% h=legend(configs);
% set(gca,'fontsize',14)
% set(h,'fontsize',14)
% ylabel('R_{d6}','fontsize',14)
