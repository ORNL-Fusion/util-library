clearvars;


fname2 = [];
% run_path = 'C:\Work\fortran\Poincare\';
% fname = 'poincare_output_120.out';  mytitle = '\phi = 120 (t=3750)'; plot_ts = 1; plot_ece=0;
% fname = 'poincare_output_120_n50.out';  mytitle = '\phi = 120 (t=3750)'; plot_ts = 1; plot_ece=0;
% fname = 'poincare_output_300.out';  mytitle = '\phi = 120 (t=4250 eqv)'; plot_ts = 1; plot_ece=0;
% fname = 'poincare_output_80.out';  mytitle = '\phi = 80 (t=3750)'; plot_ts = 0; plot_ece=1;
% fname = 'poincare_output_260.out';  mytitle = '\phi = 80  (t=4250 eqv)'; plot_ts = 0; plot_ece=1;
% fname = 'poincare_output_240.out';  mytitle = '\phi = -120 (t=3750)'; plot_ts = 1; plot_ece=0;
% fname = 'poincare_output_60.out';  mytitle = '\phi = -120 (t=4250 eqv)'; plot_ts = 1; plot_ece=0;
% fname = 'poincare_output_60_n25.out';  mytitle = '\phi = -120 (t=4250 eqv)'; plot_ts = 1; plot_ece=0; fname2 = 'psiN_min_output_60_n25.out';
% gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';

% run_path = 'C:\Work\fortran\Poincare\148712_m3dc1_t1\'; fname = 'poincare_output.out';  mytitle = '\phi = 120'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g148712.04101';
% run_path = 'C:\Work\fortran\Poincare\148712_m3dc1_t1\'; fname = 'poincare_output_240_along_TS.out';  mytitle = '\phi = 120'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g148712.04101';
% run_path = 'C:\Work\fortran\Poincare\148712_m3dc1_t1\'; fname = 'poincare_output_240_more_R.out';  mytitle = '\phi = 120'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g148712.04101'; fname2 = 'psiN_min_output_240_more_R.out';
% % run_path = 'C:\Work\fortran\Poincare\148712_m3dc1_t1\'; fname = 'poincare_output_60.out';  mytitle = '\phi = 300'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g148712.04101';
% run_path = 'C:\Work\fortran\Poincare\148712_m3dc1_t1\'; fname = 'poincare_output_60_more_R.out';  mytitle = '\phi = 300'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g148712.04101';
% run_path = 'C:\Work\fortran\Poincare\148712_m3dc1_t1\'; fname = 'poincare_output.out';  mytitle = '\phi = 120'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g148712.04101';

g=readg_g3d(gfile_name);
psiN_g = (g.psirz-g.ssimag)/(g.ssibry-g.ssimag);


fid = fopen([run_path,fname]);
dat = sscanf(fgets(fid),'%f %i %i',3);
phi_plot_deg = dat(1);
numlines = dat(2);
% % % num_pts = dat(3) - 1 ;
num_pts = dat(3)  ;

rline = zeros(numlines,num_pts);
zline = zeros(numlines,num_pts);
psiNline = zeros(numlines,num_pts);
iline = zeros(numlines,1);
for i = 1:numlines
    iline(i) = fscanf(fid,'%i',1);
    rline(i,:) = fscanf(fid,'%f',num_pts);
    zline(i,:) = fscanf(fid,'%f',num_pts);
    psiNline(i,:) = calc_psiN(g,rline(i,:),zline(i,:));
end
fclose all;

rline(isnan(psiNline)) = [];
zline(isnan(psiNline)) = [];
psiNline(isnan(psiNline)) = [];



if 0

figure; hold on; box on;
plot(rline,zline,'k.')
xlabel('R (m)','fontsize',12)
ylabel('Z (m)','fontsize',12)
set(gca,'fontsize',12)
% contour(g.r,g.z,psiN_g.',[0.5,0.6,0.7,0.8,0.9,1.0],'linewidth',2);
end
% asdf

% PLOT THETA/PSIN

tline = atan2(zline-g.zmaxis,rline-g.rmaxis);


figure; hold on; box on;
plot(tline.'/pi,psiNline.','k.')
ylabel('\psi_N','fontsize',12)
xlabel('\theta (\pi rad.)','fontsize',12)
set(gca,'fontsize',12)
title(mytitle)
% contour(g.r,g.z,psiN_g.',[0.5,0.6,0.7,0.8,0.9,1.0],'linewidth',2);
