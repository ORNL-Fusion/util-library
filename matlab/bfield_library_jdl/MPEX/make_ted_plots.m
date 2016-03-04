clearvars;

bfield_method = 1;  % 1 is Biot-Savart, 2 is exact

geo = get_Proto_geometry(0,0,0);


current_A = 3300;
% current_A = 6400;
current_B = 0;
% config = 'flat';
% helicon_current = 210;

% current_A = 6607;
% current_B = 6607;
% config = 'focus';
% helicon_current = -80;



verbose = 0;
num_test = 10;
helicon_currents = linspace(0,600,num_test);

for iconfig = 1:3
    iconfig
    if iconfig == 1
        config = 'standard';
    elseif iconfig == 2
        config = 'focus';
    elseif iconfig == 3
        config = 'flat';
    end
    for ih = 1:num_test
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
    end
end

for iconfig = 1:3
    if iconfig == 1
        config = 'standard';
    elseif iconfig == 2
        config = 'focus';
    elseif iconfig == 3
        config = 'flat';
    end
    figure; hold on; box on;
    plot(helicon_currents,1./Ru1s(iconfig,:),'b','linewidth',2)
    plot(helicon_currents,1./Ru2s(iconfig,:),'m','linewidth',2)
    plot(helicon_currents,1./Rd5s(iconfig,:),'c','linewidth',2)
    plot(helicon_currents,1./Rd6s(iconfig,:),'r','linewidth',2)
    xlabel('I_{helicon} [Amps]','fontsize',14)
    h=legend('R_{u1}','R_{u2}','R_{d5}','R_{d6}');
    ylabel('Mirror ratios','fontsize',14)
    set(gca,'fontsize',14)
    set(h,'fontsize',14)
    title(config,'fontsize',14)
end

figure; hold on; box on;
plot(helicon_currents,1./Ru1s(1,:),'b','linewidth',2)
plot(helicon_currents,1./Ru1s(2,:),'r','linewidth',2)
plot(helicon_currents,1./Ru1s(3,:),'k','linewidth',2)
xlabel('I_{helicon} [Amps]','fontsize',14)
h=legend('standard','focus','flat');
set(gca,'fontsize',14)
set(h,'fontsize',14)
ylabel('R_{u1}','fontsize',14)

figure; hold on; box on;
plot(helicon_currents,1./Rd6s(1,:),'b','linewidth',2)
plot(helicon_currents,1./Rd6s(2,:),'r','linewidth',2)
plot(helicon_currents,1./Rd6s(3,:),'k','linewidth',2)
xlabel('I_{helicon} [Amps]','fontsize',14)
h=legend('standard','focus','flat');
set(gca,'fontsize',14)
set(h,'fontsize',14)
ylabel('R_{d6}','fontsize',14)

for iconfig = 1:3
    if iconfig == 1
        config = 'standard';
    elseif iconfig == 2
        config = 'focus';
    elseif iconfig == 3
        config = 'flat';
    end
    figure; hold on; box on;
    plot(helicon_currents,r_skimmer(iconfig,:),'k','linewidth',2)
    plot(helicon_currents,r_d(iconfig,:),'r','linewidth',2)
    plot(helicon_currents,r_m(iconfig,:),'g','linewidth',2)
    plot(helicon_currents,r_u(iconfig,:),'b','linewidth',2)    

    xlabel('I_{helicon} [Amps]','fontsize',14)
    h=legend('r_{skimmer}','r_d','r_m','r_u');
    ylabel('r_skimmer','fontsize',14)
    set(gca,'fontsize',14)
    set(h,'fontsize',14)
    title(['Target 7.5, ',config],'fontsize',14)
end
