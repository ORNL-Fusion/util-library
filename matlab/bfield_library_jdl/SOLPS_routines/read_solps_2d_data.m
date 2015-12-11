clearvars;
% run_path = 'C:\Work\SOLPS\SOLPS_5_RUNS\160884\3014\test_Donly_6MW_2e19\'; fname = 'prad_data.mat'; gfile_name = 'C:\Work\M3DC1\160884\3000\g160884.03014_251';
run_path = 'C:\Work\C2_SOLPS_benchmark\152845\C2_benchmark_case1\'; fname = 'solps_data.mat'; gfile_name = 'C:\Work\C2_SOLPS_benchmark\152845\EXP_DATA\g152845.02500_516';



g = readg_g3d(gfile_name);


load([run_path,fname]);
d = solps_data.te; mytitle = 'T_e [eV]';
% d = solps_data.nedat/1e19; mytitle = 'n_e [10^{19} m^-^3]';

for i=1:solps_data.ny*solps_data.nx
   pr(:,i) = [solps_data.crx_bl(i),solps_data.crx_br(i),solps_data.crx_tr(i),solps_data.crx_tl(i)];
   pz(:,i) = [solps_data.cry_bl(i),solps_data.cry_br(i),solps_data.cry_tr(i),solps_data.cry_tl(i)];
   prcen(i) = mean(pr(:,i));
   pzcen(i) = mean(pz(:,i));
end

if 1
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

% jxa = solps_data.nx - 54 + 1;
jxa = 56;
r2d = reshape(prcen,solps_data.nx,solps_data.ny);
z2d = reshape(pzcen,solps_data.nx,solps_data.ny);
te2d = reshape(solps_data.te,solps_data.nx,solps_data.ny);
ne2d = reshape(solps_data.nedat,solps_data.nx,solps_data.ny);

% r2d = fliplr(r2d);
% te2d = fliplr(te2d);

z1d = z2d(jxa,:);
r1d = r2d(jxa,:);
te1d = te2d(jxa,:);
ne1d = ne2d(jxa,:);

plot(r1d,z1d,'k.');


figure; hold on; box on;
% contourf(r2d,z2d,te2d)
plot(r1d,te1d/1e3,'r--','linewidth',2)
xlabel('R_{OMP} [m]','fontsize',12)
ylabel('T_e [keV]','fontsize',12)
set(gca,'fontsize',12)

figure; hold on; box on;
plot(r1d,ne1d/1e19,'r--','linewidth',2)
xlabel('R_{OMP} [m]','fontsize',12)
ylabel('n_e [10^{19} m^{-3}]','fontsize',12)
set(gca,'fontsize',12)
