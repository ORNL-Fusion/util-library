clearvars;

shot = 7277;

[helicon_current,current_A,current_B,config,skimmer] = get_Proto_current(shot);
[coil,current] = build_Proto_coils(helicon_current,current_A,current_B,config);
% [rr_cm_IR,dd_cm_IR,plasma_radius_cm] = plot_IR_data_raw(shot,1,0,-2.5);
geo = get_Proto_geometry(0,0,skimmer);

bfield.coil = coil;
bfield.current = current;
bfield.type = 'just_coils';

num_lines = 20;

tic
% Forward and reverse lines from target
rr = linspace(1e-3,0.03,num_lines);
zz = geo.target.z*ones(size(rr));
L = geo.target.z - 0.5;
dl = -0.01; nsteps = round(abs(L/dl)); dl = sign(dl)*L/nsteps;
phistart = zeros(size(rr));
f2a = follow_fieldlines_rzphi_dz(bfield,rr,zz(1),phistart,dl,nsteps);
f2a = clip_fl_at_vessel(f2a,geo);
% lines from helicon edges
% back
L = geo.target.z - geo.helicon.z1;
dl = 0.01; nsteps = round(abs(L/dl)); dl = sign(dl)*L/nsteps;
f3a = follow_fieldlines_rzphi_dz(bfield,geo.helicon.r,geo.helicon.z1,0,dl,nsteps);
f3a = clip_fl_at_vessel(f3a,geo);
L = geo.helicon.z1 - 0.5;
dl = -0.01; nsteps = round(abs(L/dl)); dl = sign(dl)*L/nsteps;
f3b = follow_fieldlines_rzphi_dz(bfield,geo.helicon.r,geo.helicon.z1,0,dl,nsteps);
f3b = clip_fl_at_vessel(f3b,geo);
%front
L = geo.target.z - geo.helicon.z2;
dl = 0.01; nsteps = round(abs(L/dl)); dl = sign(dl)*L/nsteps;
f4a = follow_fieldlines_rzphi_dz(bfield,geo.helicon.r,geo.helicon.z2,0,dl,nsteps);
f4a = clip_fl_at_vessel(f4a,geo);
L = geo.helicon.z2 - 0.5;
dl = -0.01; nsteps = round(abs(L/dl)); dl = sign(dl)*L/nsteps;
f4b = follow_fieldlines_rzphi_dz(bfield,geo.helicon.r,geo.helicon.z2,0,dl,nsteps);
f4b = clip_fl_at_vessel(f4b,geo);
%skimmer
if skimmer
    L = geo.target.z - geo.skimmer.zmid;
    dl = 0.01; nsteps = round(abs(L/dl)); dl = sign(dl)*L/nsteps;
    f5a = follow_fieldlines_rzphi_dz(bfield,geo.skimmer.r,geo.skimmer.zmid,0,dl,nsteps);
    f5a = clip_fl_at_vessel(f5a,geo);
    L = geo.skimmer.zmid - 0.5;
    dl = -0.01; nsteps = round(abs(L/dl)); dl = sign(dl)*L/nsteps;
    f5b = follow_fieldlines_rzphi_dz(bfield,geo.skimmer.r,geo.skimmer.zmid,0,dl,nsteps);
    f5b = clip_fl_at_vessel(f5b,geo);
end
toc


ilcfs = find(~isnan(f2a.z(end,:)),1,'last');


figure; hold on; box on;
plot(f2a.z,f2a.r,'b','linewidth',2); 
plot(f2a.z(:,ilcfs),f2a.r(:,ilcfs),'r')
plot(f3a.z,f3a.r,'r','linewidth',2); plot(f3b.z,f3b.r,'r','linewidth',2)
plot(f4a.z,f4a.r,'r','linewidth',2); plot(f4b.z,f4b.r,'r','linewidth',2)
plot(f5a.z,f5a.r,'c','linewidth',2); plot(f5b.z,f5b.r,'c','linewidth',2)
set(gca,'fontsize',14)
xlabel('Z [m]','fontsize',14)
ylabel('R [m]','fontsize',14)
title(['Shot ',num2str(shot)])
get_Proto_geometry(1,0,skimmer);
axis([0.5,3.5,0,0.2])


% plot(

