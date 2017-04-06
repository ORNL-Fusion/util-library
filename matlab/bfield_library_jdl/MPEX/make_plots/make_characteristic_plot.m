function make_characteristic_plot(shot)

% shot = 7477;

ZMIN = 1;

[helicon_current,current_A,current_B,config,skimmer] = get_Proto_current(shot);
% helicon_current = 100;
% current_A = 3300;
% current_B = 0;
% config = 'standard';
% skimmer = 1;
% 
% shot_tmp{1} = helicon_current;
% shot_tmp{2} = current_A;
% shot_tmp{3} = current_B;
% shot_tmp{4} = config;
% shot_tmp{5} = skimmer;
% f = find_lcfs(shot_tmp);

[coil,current] = build_Proto_coils(helicon_current,current_A,current_B,config);
% [rr_cm_IR,dd_cm_IR,plasma_radius_cm] = plot_IR_data_raw(shot,1,0,-2.5);
geo = get_Proto_geometry(0,0,skimmer);
f_lcfs = find_lcfs(shot);
% f_lcfs = find_lcfs(shot_tmp,0,6,0.01);

bfield.coil = coil;
bfield.current = current;
bfield.type = 'just_coils';

num_lines = 10;

tic
% Forward and reverse lines from helicon center
rr = linspace(1e-3,geo.helicon.r,num_lines);
zz = geo.helicon.zmid*ones(size(rr));
L = geo.target.z - geo.helicon.zmid;
dl = 0.01; nsteps = round(abs(L/dl)); dl = sign(dl)*L/nsteps;
phistart = zeros(size(rr));
f2a = follow_fieldlines_rzphi_dz(bfield,rr,zz(1),phistart,dl,nsteps);
f2a = clip_fl_at_vessel(f2a,geo);
L = geo.helicon.zmid - ZMIN;
dl = -0.01;
nsteps = round(abs(L/dl)); dl = sign(dl)*L/nsteps;
f2b = follow_fieldlines_rzphi_dz(bfield,rr,zz(1),phistart,dl,nsteps);
f2b = clip_fl_at_vessel(f2b,geo);
% lines from helicon edges
% back
L = geo.target.z - geo.helicon.z1;
dl = 0.01; nsteps = round(abs(L/dl)); dl = sign(dl)*L/nsteps;
f3a = follow_fieldlines_rzphi_dz(bfield,geo.helicon.r,geo.helicon.z1,0,dl,nsteps);
f3a = clip_fl_at_vessel(f3a,geo);
L = geo.helicon.z1 - ZMIN;
dl = -0.01; nsteps = round(abs(L/dl)); dl = sign(dl)*L/nsteps;
f3b = follow_fieldlines_rzphi_dz(bfield,geo.helicon.r,geo.helicon.z1,0,dl,nsteps);
f3b = clip_fl_at_vessel(f3b,geo);
%front
L = geo.target.z - geo.helicon.z2;
dl = 0.01; nsteps = round(abs(L/dl)); dl = sign(dl)*L/nsteps;
f4a = follow_fieldlines_rzphi_dz(bfield,geo.helicon.r,geo.helicon.z2,0,dl,nsteps);
f4a = clip_fl_at_vessel(f4a,geo);
L = geo.helicon.z2 - ZMIN;
dl = -0.01; nsteps = round(abs(L/dl)); dl = sign(dl)*L/nsteps;
f4b = follow_fieldlines_rzphi_dz(bfield,geo.helicon.r,geo.helicon.z2,0,dl,nsteps);
f4b = clip_fl_at_vessel(f4b,geo);
%skimmer
if skimmer
    L = geo.target.z - geo.skimmer.zmid;
    dl = 0.01; nsteps = round(abs(L/dl)); dl = sign(dl)*L/nsteps;
    f5a = follow_fieldlines_rzphi_dz(bfield,geo.skimmer.r,geo.skimmer.zmid,0,dl,nsteps);
    f5a = clip_fl_at_vessel(f5a,geo);
    L = geo.skimmer.zmid - ZMIN;
    dl = -0.01; nsteps = round(abs(L/dl)); dl = sign(dl)*L/nsteps;
    f5b = follow_fieldlines_rzphi_dz(bfield,geo.skimmer.r,geo.skimmer.zmid,0,dl,nsteps);
    f5b = clip_fl_at_vessel(f5b,geo);
end
toc


figure; hold on; box on;
plot(f2a.z,f2a.r,'b','linewidth',2); plot(f2b.z,f2b.r,'b','linewidth',2)
plot(f_lcfs.z,f_lcfs.r,'r','linewidth',2)
plot(f3a.z,f3a.r,'c','linewidth',2); plot(f3b.z,f3b.r,'c','linewidth',2)
plot(f4a.z,f4a.r,'c','linewidth',2); plot(f4b.z,f4b.r,'c','linewidth',2)
if skimmer
    plot(f5a.z,f5a.r,'m','linewidth',2); plot(f5b.z,f5b.r,'m','linewidth',2)
end
set(gca,'fontsize',14)
xlabel('Z [m]','fontsize',14)
ylabel('R [m]','fontsize',14)
title(sprintf('Shot %d, I_H=%3.0f A, I_A=%4.0f A',shot,helicon_current,current_A))
get_Proto_geometry(1,0,skimmer);
axis([0.5,3.5,0,0.15])


