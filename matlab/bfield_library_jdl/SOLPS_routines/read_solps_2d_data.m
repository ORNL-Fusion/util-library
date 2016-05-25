function solps_data = read_solps_2d_data(run_path,fname,g,just_return_data)
if nargin == 0
% run_path = 'C:\Work\SOLPS\SOLPS_5_RUNS\160884\3014\test_Donly_6MW_2e19\'; fname = 'prad_data.mat'; gfile_name = 'C:\Work\M3DC1\160884\3000\g160884.03014_251';
% run_path = 'C:\Work\C2_SOLPS_benchmark\152845\C2_benchmark_case1\'; fname = 'solps_data.mat'; gfile_name = 'C:\Work\C2_SOLPS_benchmark\152845\EXP_DATA\g152845.02500_516';

% run_path = 'C:\Work\CMOD\PSI2016\1150625014\1009_D+Ne\P2MW_n_1e20_puff_3.3e19'; fname = 'data_2d.mat'; gfile_name = 'C:\Work\CMOD\PSI2016\1150625014\1009_D+Ne\baserun\g1150625014.01009_983';
run_path = 'C:\Work\CMOD\PSI2016\1150625014\SOLPS5\1009_D+N\P2MW_n_2e20_sxp_no_puff_leakage_moresxp_adjust_1.2\'; fname = 'data_2d.mat'; gfile_name = 'C:\Work\CMOD\PSI2016\1150625014\1009_D+Ne\baserun\g1150625014.01009_983';
just_return_data = 0;
g = readg_g3d(gfile_name);
end
jpol = 74;
% jpol = 86
% jpol = 95

plot_geo_pts = 1;

% 

plot_profs = 1; rad_type = 1;  % 1 is psiN, 0 is Rmid

dat = load(fullfile(run_path,fname));
solps_data = dat.solps_data;

for i=1:solps_data.ny*solps_data.nx
   pr(:,i) = [solps_data.crx_bl(i),solps_data.crx_br(i),solps_data.crx_tr(i),solps_data.crx_tl(i)];
   pz(:,i) = [solps_data.cry_bl(i),solps_data.cry_br(i),solps_data.cry_tr(i),solps_data.cry_tl(i)];
   prcen(i) = mean(pr(:,i));
   pzcen(i) = mean(pz(:,i));
   prleft(i) = mean(pr([1,4],i));
   pzleft(i) = mean(pz([1,4],i));
   prright(i) = mean(pr([2,3],i));
   pzright(i) = mean(pz([2,3],i));
end

r2d_left = reshape(prleft,solps_data.nx,solps_data.ny);
z2d_left = reshape(pzleft,solps_data.nx,solps_data.ny);


solps_data.r2d_left = r2d_left;
solps_data.z2d_left = z2d_left;
solps_data.pr = pr;
solps_data.pz = pz;

if just_return_data
    return;
end
fprintf('Using jpol = %d\n',jpol)

d = solps_data.te; mytitle = 'T_e [eV]';
% d = solps_data.nedat/1e19; mytitle = 'n_e [10^{19} m^-^3]';
% d = (solps_data.dens_n0 + solps_data.dens_n1 +solps_data.dens_n2 +solps_data.dens_n3+solps_data.dens_n4+solps_data.dens_n5+solps_data.dens_n6+solps_data.dens_n7); mytitle = 'n N total [m^{-3}]'
% d = solps_data.prad_dfluid/1e6; mytitle = 'prad D fluid[MW/m^3]';
% d = solps_data.prad_dneut/1e6; mytitle = 'prad D neutral [MW/m^3]';
% d = solps_data.prad_dbrem/1e6; mytitle = 'prad D brem [MW/m^3]';
% d = (solps_data.prad_dfluid + solps_data.prad_dneut + solps_data.prad_dbrem)/1e6; mytitle = 'prad D tot [MW/m^3]';

% d = solps_data.prad_nefluid/1e6; mytitle = 'prad Ne fluid[MW/m^3]';
% d = solps_data.prad_neneut/1e6; mytitle = 'prad Ne neut[MW/m^3]';
% d = (solps_data.prad_nefluid + solps_data.prad_neneut + solps_data.prad_nebrem)/1e6; mytitle = 'prad Ne tot [MW/m^3]';
% d = solps_data.vol; mytitle = 'Cell volume [m^{-3}]';




if plot_geo_pts
    figure; hold on; box on;
    data_plot = d;
    patch(pr,pz,data_plot,'edgecolor','none')
    colorbar;
    plot(g.lim(1,g.lim(1,:)>0),g.lim(2,g.lim(1,:)>0),'k')
    axis tight;
    colorbar('fontsize',14)
    xlabel('R [m]','fontsize',14)
    ylabel('Z [m]','fontsize',14)
    set(gca,'fontsize',14)
    title(mytitle);
end

% jpol = solps_data.nx - 54 + 1;
% jpol = 56;


te2d = reshape(solps_data.te,solps_data.nx,solps_data.ny);
ne2d = reshape(solps_data.nedat,solps_data.nx,solps_data.ny);
r2d = reshape(prcen,solps_data.nx,solps_data.ny);
z2d = reshape(pzcen,solps_data.nx,solps_data.ny);
r2d_right = reshape(prright,solps_data.nx,solps_data.ny);
z2d_right = reshape(pzright,solps_data.nx,solps_data.ny);

% r2d = fliplr(r2d);
% te2d = fliplr(te2d);

z1d = z2d(jpol,:);
r1d = r2d(jpol,:);
z1d_left = z2d_left(jpol,:);
r1d_left = r2d_left(jpol,:);
z1d_right = z2d_right(jpol,:);
r1d_right = r2d_right(jpol,:);

te1d = te2d(jpol,:);
ne1d = ne2d(jpol,:);

plot(r1d,z1d,'k.');
plot(r1d_left,z1d_left,'co');
plot(r1d_right,z1d_right,'mo');


psiN_rad = calc_psiN(g,r1d_left,z1d_left);
psiN_cen = calc_psiN(g,r1d,z1d);
psiN_left = calc_psiN(g,r1d_left,z1d_left);
psiN_right = calc_psiN(g,r1d_right,z1d_right);
% figure; hold on;
% plot(psiN_cen)
% plot(psiN_right)
% plot(psiN_left)


if plot_profs
    figure; hold on; box on;
    if rad_type == 0
        plot(r1d,te1d,'r--','linewidth',2)
        xlabel('R_{OMP} [m]','fontsize',12)
    elseif rad_type == 1        
        plot(psiN_rad,te1d,'r--','linewidth',2)
        xlabel('\Psi_N','fontsize',12)
    else
        error('Unknown rad_type')
    end
    
    ylabel('T_e [eV]','fontsize',12)
    set(gca,'fontsize',12)
    
    figure; hold on; box on;
    if rad_type == 0
        plot(r1d,ne1d/1e20,'r--','linewidth',2)
        xlabel('R_{OMP} [m]','fontsize',12)
    elseif rad_type == 1
        plot(psiN_rad,ne1d/1e20,'r--','linewidth',2)
        xlabel('\Psi_N','fontsize',12)
    else
        error('Unknown rad_type')
    end
    
    
    ylabel('n_e [10^{20} m^{-3}]','fontsize',12)
    set(gca,'fontsize',12)
end