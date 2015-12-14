clearvars;


plot_glim = 0;

fname_psimin = []; fname_psimin2 = [];
% run_path = 'C:\Work\fortran\Poincare\';
% fname = 'poincare_output_120.out';  mytitle = '\phi = 120 (t=3750)'; plot_ts = 1; plot_ece=0;
% fname = 'poincare_output_120_n50.out';  mytitle = '\phi = 120 (t=3750)'; plot_ts = 1; plot_ece=0;
% fname = 'poincare_output_300.out';  mytitle = '\phi = 120 (t=4250 eqv)'; plot_ts = 1; plot_ece=0;
% fname = 'poincare_output_80.out';  mytitle = '\phi = 80 (t=3750)'; plot_ts = 0; plot_ece=1;
% fname = 'poincare_output_260.out';  mytitle = '\phi = 80  (t=4250 eqv)'; plot_ts = 0; plot_ece=1;
% fname = 'poincare_output_240.out';  mytitle = '\phi = -120 (t=3750)'; plot_ts = 1; plot_ece=0;
% fname = 'poincare_output_60.out';  mytitle = '\phi = -120 (t=4250 eqv)'; plot_ts = 1; plot_ece=0;
% fname = 'poincare_output_60_n25.out';  mytitle = '\phi = -120 (t=4250 eqv)'; plot_ts = 1; plot_ece=0; fname_psimin = 'psiN_min_output_60_n25.out';fname_psimin = [];
% gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';

% run_path = 'C:\Work\fortran\Poincare\148712_m3dc1_t1\'; fname = 'poincare_output_240.out';  mytitle = '\phi = 120'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g148712.04101';
% run_path = 'C:\Work\fortran\Poincare\148712_m3dc1_t1\'; fname = 'poincare_output_240_along_TS.out';  mytitle = '\phi = 120'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g148712.04101';
% run_path = 'C:\Work\fortran\Poincare\148712_m3dc1_t1\'; fname = 'poincare_output_240_more_R.out';  mytitle = '\phi = 120'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g148712.04101'; fname_psimin = 'psiN_min_output_240_more_R.out';
% run_path = 'C:\Work\fortran\Poincare\148712_m3dc1_t1\'; fname = 'poincare_output_60.out';  mytitle = '\phi = 300'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g148712.04101';
% run_path = 'C:\Work\fortran\Poincare\148712_m3dc1_t1\'; fname = 'poincare_output_60_more_R.out';  mytitle = '\phi = 300'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g148712.04101';
% run_path = 'C:\Work\fortran\Poincare\148712_m3dc1_t1\'; fname = 'poincare_output.out';  mytitle = '\phi = 120'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g148712.04101';

% run_path = 'C:\Work\fortran\Poincare\154929\rmp\';fname = 'poincare_output_60_n500.out';  fname2 = 'poincare_output2_60_n500.out';mytitle = '\phi = 60 (t=4250 eqv)'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output_60_n500.out';fname_psimin2 = 'psiN_min_output2_60_n500.out';shot = 154921;times = 4250;
% run_path = 'C:\Work\fortran\Poincare\154929\rmp\';fname = 'poincare_output_60_n500_fake_diag.out';  fname2 = 'poincare_output2_60_n500_fake_diag.out';mytitle = '\phi = 60 (t=4250 eqv)'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output_60_n500_fake_diag.out';fname_psimin2 = 'psiN_min_output2_60_n500_fake_diag.out';shot = 154921;times = 4250;
% run_path = 'C:\Work\fortran\Poincare\154929\rmp\';fname = 'poincare_output__60_n100_fake_diag2.out';  fname2 = 'poincare_output2__60_n100_fake_diag2.out';mytitle = '\phi = 60 (t=4250 eqv)'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output__60_n100_fake_diag2.out';fname_psimin2 = 'psiN_min_output2__60_n100_fake_diag2.out';shot = 154921;times = 4250;
% run_path = 'C:\Work\fortran\Poincare\154929\rmp\';fname = 'poincare_output_240_n500.out';  fname2 = 'poincare_output2_240_n500.out';mytitle = '\phi = -120 (t=3750)'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output_240_n500.out';fname_psimin2 = 'psiN_min_output2_240_n500.out';shot = 154921;times = 3750;
% run_path = 'C:\Work\fortran\Poincare\154929\rmp\';fname = 'poincare_output_100_n500.out';  fname2 = 'poincare_output2_100_n500.out';mytitle = '\phi = 280 (t=4250 eqv)'; plot_ts = 0; plot_ece=1;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output_100_n500.out';fname_psimin2 = 'psiN_min_output2_100_n500.out';shot = 154921;times = 4250;
% run_path = 'C:\Work\fortran\Poincare\154929\rmp\';fname = 'poincare_output_280_n500.out';  fname2 = 'poincare_output2_280_n500.out';mytitle = '\phi = -80 (t=3750)'; plot_ts = 0; plot_ece=1;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output_280_n500.out';fname_psimin2 = 'psiN_min_output2_280_n500.out';shot = 154921;times = 3750;

% run_path = 'C:\Work\fortran\Poincare\154929/m3dc1_t1\';fname = 'poincare_output_240_n500.out';  fname2 = 'poincare_output2_240_n500.out';mytitle = '\phi = -120 (t=3750)'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output_240_n500.out';fname_psimin2 = 'psiN_min_output2_240_n500.out';shot = 154921;times = 3750;
% run_path = 'C:\Work\fortran\Poincare\154929_m3dc1_t1\';fname = 'poincare_output_60_n500.out';  fname2 = 'poincare_output2_60_n500.out';mytitle = '\phi = 60 (t=4250 eqv)'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output_60_n500.out';fname_psimin2 = 'psiN_min_output2_60_n500.out';shot = 154921;times = 4250;
% run_path = 'C:\Work\fortran\Poincare\154929\m3dc1_t1\';fname = 'poincare_output_100_n500.out';  fname2 = 'poincare_output2_100_n500.out';mytitle = '\phi = 280 (t=4250 eqv)'; plot_ts = 0; plot_ece=1;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output_100_n500.out';fname_psimin2 = 'psiN_min_output2_100_n500.out';shot = 154921;times = 4250;

% run_path = 'C:\Work\fortran\Poincare\154929_m3dc1_t1_fb\';fname = 'poincare_output_60_n500.out';  fname2 = 'poincare_output2_60_n500.out';mytitle = '\phi = 60 (t=4250 eqv)'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output_60_n500.out';fname_psimin2 = 'psiN_min_output2_60_n500.out';shot = 154921;times = 4250;
% run_path = 'C:\Work\fortran\Poincare\154929_m3dc1_t1_fb\';fname = 'poincare_output_240_n500.out';  fname2 = 'poincare_output2_240_n500.out';mytitle = '\phi = -120 (t=3750)'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output_240_n500.out';fname_psimin2 = 'psiN_min_output2_240_n500.out';shot = 154921;times = 3750;

% run_path = 'C:\Work\fortran\Poincare\154929\';fname = 'poincare_output_240AS.out';  fname2 = 'poincare_output2_240AS.out';mytitle = '\phi = -120 (t=3750)'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output_240AS.out';fname_psimin2 = 'psiN_min_output2_240AS.out';

% run_path = 'C:\Work\fortran\Poincare\154929\';fname = 'poincare_output_240_n250.out';  fname2 = 'poincare_output2_240_n250.out';mytitle = '\phi = -120 (t=3750)'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output_240_n250.out';fname_psimin2 = 'psiN_min_output2_240_n250.out';shot = 154921;times = 3750;
% run_path = 'C:\Work\fortran\Poincare\154929_current_scan\';fname = 'poincare_output_240_s1.25.out';  fname2 = 'poincare_output2_240_s1.25.out';mytitle = '\phi = -120 (t=3750)'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output_240_s1.25.out';fname_psimin2 = 'psiN_min_output2_240_s1.25.out';shot = 154921;times = 3750;
% run_path = 'C:\Work\fortran\Poincare\154929_current_scan\';fname = 'poincare_output_240_s2.out';  fname2 = 'poincare_output2_240_s2.out';mytitle = '\phi = -120 (t=3750)'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output_240_s2.out';fname_psimin2 = 'psiN_min_output2_240_s2.out';shot = 154921;times = 3750;
% run_path = 'C:\Work\fortran\Poincare\154929_current_scan\';fname = 'poincare_output_240_s4.out';  fname2 = 'poincare_output2_240_s4.out';mytitle = '\phi = -120 (t=3750)'; plot_ts = 1; plot_ece=0;gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';fname_psimin = 'psiN_min_output_240_s4.out';fname_psimin2 = 'psiN_min_output2_240_s4.out';shot = 154921;times = 3750;

% run_path = 'C:\Work\fortran\Poincare\160884\04015\rmp\';suffix = 'dts';fname = ['poincare_output_',suffix,'.out'];  fname2 = ['poincare_output2_',suffix,'.out'];mytitle = '\phi = 0'; plot_ts = 0; plot_ece=0;gfile_name = 'C:\Work\DIII-D\160884\g160884.04015_915';fname_psimin = ['psiN_min_output_',suffix,'.out'];fname_psimin2 = ['psiN_min_output2_',suffix,'.out'];
% run_path = 'C:\Work\fortran\Poincare\160884\04015\rmp\';suffix = 'Z0_fine';fname = ['poincare_output_',suffix,'.out'];  fname2 = ['poincare_output2_',suffix,'.out'];mytitle = '\phi = 0'; plot_ts = 0; plot_ece=0;gfile_name = 'C:\Work\DIII-D\160884\g160884.04015_915';fname_psimin = ['psiN_min_output_',suffix,'.out'];fname_psimin2 = ['psiN_min_output2_',suffix,'.out'];

% run_path = 'C:\Work\fortran\Poincare\160884\03014\m3dc1_full\'; suffix = '0deg_n100_moreR';fname = ['poincare_output_',suffix,'.out'];  fname2 = ['poincare_output2_',suffix,'.out'];mytitle = '\phi = 0'; plot_ts = 0; plot_ece=0;gfile_name = 'C:\Work\M3DC1\160884\3000\g160884.03014_251';fname_psimin = ['psiN_min_output_',suffix,'.out'];fname_psimin2 = ['psiN_min_output2_',suffix,'.out'];
% run_path = 'C:\Work\fortran\Poincare\160884\03014\m3dc1_g+pert\'; suffix = '0deg_n100';fname = ['poincare_output_',suffix,'.out'];  fname2 = ['poincare_output2_',suffix,'.out'];mytitle = '\phi = 0'; plot_ts = 0; plot_ece=0;gfile_name = 'C:\Work\M3DC1\160884\3000\g160884.03014_251';fname_psimin = ['psiN_min_output_',suffix,'.out'];fname_psimin2 = ['psiN_min_output2_',suffix,'.out'];
% run_path = 'C:\Work\fortran\Poincare\160884\03014\g3d_as\'; suffix = '0deg_n100';fname = ['poincare_output_',suffix,'.out'];  fname2 = ['poincare_output2_',suffix,'.out'];mytitle = '\phi = 0'; plot_ts = 0; plot_ece=0;gfile_name = 'C:\Work\M3DC1\160884\3000\g160884.03014_251';fname_psimin = ['psiN_min_output_',suffix,'.out'];fname_psimin2 = ['psiN_min_output2_',suffix,'.out'];
% run_path = 'C:\Work\fortran\Poincare\160884\03014\m3dc1_as\'; suffix = '0deg_n100';fname = ['poincare_output_',suffix,'.out'];  fname2 = ['poincare_output2_',suffix,'.out'];mytitle = '\phi = 0'; plot_ts = 0; plot_ece=0;gfile_name = 'C:\Work\M3DC1\160884\3000\g160884.03014_251';fname_psimin = ['psiN_min_output_',suffix,'.out'];fname_psimin2 = ['psiN_min_output2_',suffix,'.out'];
% run_path = 'C:\Work\fortran\Poincare\160884\03014\m3dc1_g+pert_t0\'; suffix = '0deg_n100';fname = ['poincare_output_',suffix,'.out'];  fname2 = ['poincare_output2_',suffix,'.out'];mytitle = '\phi = 0'; plot_ts = 0; plot_ece=0;gfile_name = 'C:\Work\M3DC1\160884\3000\g160884.03014_251';fname_psimin = ['psiN_min_output_',suffix,'.out'];fname_psimin2 = ['psiN_min_output2_',suffix,'.out'];

% run_path = 'C:\Work\fortran\Poincare\160884\03014\m3dc1_as\'; suffix = '0deg_n100_0fac';fname = ['poincare_output_',suffix,'.out'];  fname2 = ['poincare_output2_',suffix,'.out'];mytitle = '\phi = 0'; plot_ts = 0; plot_ece=0;gfile_name = 'C:\Work\M3DC1\160884\3000\g160884.03014_251';fname_psimin = ['psiN_min_output_',suffix,'.out'];fname_psimin2 = ['psiN_min_output2_',suffix,'.out'];
% run_path = 'C:\Work\fortran\Poincare\160884\03014\m3dc1_g+pert_t0\'; suffix = '0deg_n100_0fac';fname = ['poincare_output_',suffix,'.out'];  fname2 = ['poincare_output2_',suffix,'.out'];mytitle = '\phi = 0'; plot_ts = 0; plot_ece=0;gfile_name = 'C:\Work\M3DC1\160884\3000\g160884.03014_251';fname_psimin = ['psiN_min_output_',suffix,'.out'];fname_psimin2 = ['psiN_min_output2_',suffix,'.out'];
 
% run_path = 'C:\Work\fortran\Poincare\160884\03014\m3dc1_g+pert\'; suffix = '0deg_n250';fname = ['poincare_output_',suffix,'.out'];  fname2 = ['poincare_output2_',suffix,'.out'];mytitle = '\phi = 0'; plot_ts = 0; plot_ece=0;gfile_name = 'C:\Work\M3DC1\160884\3000\g160884.03014_251';fname_psimin = ['psiN_min_output_',suffix,'.out'];fname_psimin2 = ['psiN_min_output2_',suffix,'.out'];

% run_path = 'C:\Work\fortran\Poincare\160884\5009\m3dc1_full\'; suffix = 'test';fname = ['poincare_output_',suffix,'.out'];  fname2 = ['poincare_output2_',suffix,'.out'];mytitle = '\phi = 0'; plot_ts = 0; plot_ece=0;gfile_name = 'C:\Work\M3DC1\160884\5000\g160884.05009_537';fname_psimin = ['psiN_min_output_',suffix,'.out'];fname_psimin2 = ['psiN_min_output2_',suffix,'.out'];
run_path = 'C:\Work\fortran\Poincare\160884\5009\m3dc1_g+pert\'; suffix = 'test';
fname = ['poincare_output_',suffix,'.out'];  fname2 = ['poincare_output2_',suffix,'.out'];mytitle = '\phi = 0'; plot_ts = 0; plot_ece=0;gfile_name = 'C:\Work\M3DC1\160884\5000\g160884.05009_537';fname_psimin = ['psiN_min_output_',suffix,'.out'];fname_psimin2 = ['psiN_min_output2_',suffix,'.out'];
run_path ='C:\Work\'; suffix = 'test';
psiN_max_eval = 0.98;

g=readg_g3d(gfile_name);
psiN_g = (g.psirz-g.ssimag)/(g.ssibry-g.ssimag);

max_file = 1 + ~isempty(fname2);
for ifile = 1:max_file
    if ifile == 1
        fname_tmp = fname;
    else
        fname_tmp = fname2;
    end
    fid = fopen([run_path,fname_tmp]);
    if fid == -1
        fprintf('ERROR: cannot open %s\n',[run_path,fname_tmp]);
        error('bad')
    end
    dat = sscanf(fgets(fid),'%f %i %i',3);    
    phi_plot_deg = dat(1);
    numlines = dat(2);
    num_pts = dat(3)  ;
    
    rline{ifile} = zeros(numlines,num_pts);
    zline{ifile} = zeros(numlines,num_pts);
    psiNline{ifile} = zeros(numlines,num_pts);
    iline{ifile} = zeros(numlines);
    for i = 1:numlines
        iline{ifile}(i) = fscanf(fid,'%i',1);
        rline{ifile}(i,:) = fscanf(fid,'%f',num_pts);
        zline{ifile}(i,:) = fscanf(fid,'%f',num_pts);
        psiNline{ifile}(i,:) = calc_psiN(g,rline{ifile}(i,:),zline{ifile}(i,:),1);
    end
    fclose all;
    
%     itmp = isnan(psiNline{ifile});
    
    rline{ifile}(isnan(psiNline{ifile})) = NaN;
    zline{ifile}(isnan(psiNline{ifile})) =NaN;
    psiNline{ifile}(isnan(psiNline{ifile})) = NaN;
end


%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
% PLOT RZ poincare
%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------

if 1    
    figure; hold on; box on;
    for ifile = 1:max_file
        plot(rline{ifile},zline{ifile},'k.','markersize',2)
    end
    xlabel('R (m)','fontsize',12)
    ylabel('Z (m)','fontsize',12)
    set(gca,'fontsize',12)
    % contour(g.r,g.z,psiN_g.',[0.5,0.6,0.7,0.8,0.9,1.0],'linewidth',2);
    if plot_ts
        R(1:40)=1.94;
        Z=[0.0430,0.1720,0.2810,0.3320,0.3590,0.3910,0.4210,0.4770,0.5310,0.5630,0.5790,0.5950,0.6080,0.6230,0.6350,0.6410,0.6480,0.6540,0.6610,0.6680,...
            0.6740,0.6810,0.6870,0.6930,0.7010,0.7060,0.7130,0.7180,0.7250,0.7310,0.7370,0.7430,0.7490,0.7550,0.7630,0.7760,0.7880,0.8000,0.8130,0.8270];
        plot(R,Z,'ro','linewidth',2)
    end
    if plot_ece
        Z2(1:40) = 0.004;
        R2 = [2.2399,2.2175,2.1951,2.1724,2.1496,2.1269,2.1045,2.0824,2.0608,2.0398,2.0193,1.9994,1.9801,1.9613,1.9430,1.9252,1.9252,1.9077,1.8905,1.8737,...
            1.8570,1.8406,1.8244,1.8084,1.7925,1.7768,1.7612,1.7458,1.7306,1.7154,1.7005,1.6856,1.6563,1.6276,1.5994,1.5717,1.5445,1.5178,1.4916,1.4658];
        T2 = atan2(Z2-g.zmaxis,R2-g.rmaxis);
        PN2 = (g.ip_sign*get_psi_bicub(g,R2,Z2)-g.ssimag)/(g.ssibry-g.ssimag);
        plot(R2,Z2,'bo','linewidth',3)
    end
    if plot_glim
        plot(g.lim(1,:),g.lim(2,:),'k','linewidth',2)
    end
end

%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
% PLOT THETA/PSIN
%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------

for ifile = 1:max_file
    tline{ifile} = atan2(zline{ifile}-g.zmaxis,rline{ifile}-g.rmaxis);
end


figure; hold on; box on;
for ifile = 1:max_file
    plot(psiNline{ifile}.',tline{ifile}.'/pi,'k.','markersize',2)
end
xlabel('\psi_N','fontsize',12)
ylabel('\theta (\pi rad.)','fontsize',12)
set(gca,'fontsize',12)
title(mytitle)
% contour(g.r,g.z,psiN_g.',[0.5,0.6,0.7,0.8,0.9,1.0],'linewidth',2);

if plot_ts
    R(1:40)=1.94;
    Z=[0.0430,0.1720,0.2810,0.3320,0.3590,0.3910,0.4210,0.4770,0.5310,0.5630,0.5790,0.5950,0.6080,0.6230,0.6350,0.6410,0.6480,0.6540,0.6610,0.6680,...
        0.6740,0.6810,0.6870,0.6930,0.7010,0.7060,0.7130,0.7180,0.7250,0.7310,0.7370,0.7430,0.7490,0.7550,0.7630,0.7760,0.7880,0.8000,0.8130,0.8270];
    T = atan2(Z-g.zmaxis,R-g.rmaxis);
    PN = (g.ip_sign*get_psi_bicub(g,R,Z)-g.ssimag)/(g.ssibry-g.ssimag);
    plot(PN,T/pi,'ro','linewidth',3)
end

if plot_ece
    Z2(1:40) = 0.004;
    R2 = [2.2399,2.2175,2.1951,2.1724,2.1496,2.1269,2.1045,2.0824,2.0608,2.0398,2.0193,1.9994,1.9801,1.9613,1.9430,1.9252,1.9252,1.9077,1.8905,1.8737,...
        1.8570,1.8406,1.8244,1.8084,1.7925,1.7768,1.7612,1.7458,1.7306,1.7154,1.7005,1.6856,1.6563,1.6276,1.5994,1.5717,1.5445,1.5178,1.4916,1.4658];
    T2 = atan2(Z2-g.zmaxis,R2-g.rmaxis);
    PN2 = (g.ip_sign*get_psi_bicub(g,R2,Z2)-g.ssimag)/(g.ssibry-g.ssimag);
    plot(PN2,T2/pi,'bo','linewidth',3)
end

if plot_glim
    rlim = g.lim(1,:);
    zlim = g.lim(2,:);
    tlim = atan2(zlim-g.zmaxis,rlim-g.rmaxis);
    psilim = (g.ip_sign*get_psi_bicub(g,rlim,zlim)-g.ssimag)/(g.ssibry-g.ssimag);    
%     plot(psilim,tlim/pi,'c-','linewidth',3)
    
    ninterp_lim = 100;
    for i = 1:size(g.lim,2)-1
        rinterp_lim = linspace(g.lim(1,i),g.lim(1,i+1),ninterp_lim);
        zinterp_lim = linspace(g.lim(2,i),g.lim(2,i+1),ninterp_lim);
        psilim_interp = (g.ip_sign*get_psi_bicub(g,rinterp_lim,zinterp_lim)-g.ssimag)/(g.ssibry-g.ssimag);    
        tlim_interp = atan2(zinterp_lim-g.zmaxis,rinterp_lim-g.rmaxis);
        plot(psilim_interp,tlim_interp/pi,'bx','linewidth',3)
    end
        rinterp_lim = linspace(g.lim(1,end),g.lim(1,end),ninterp_lim);
        zinterp_lim = linspace(g.lim(2,1),g.lim(2,1),ninterp_lim);
        psilim_interp = (g.ip_sign*get_psi_bicub(g,rinterp_lim,zinterp_lim)-g.ssimag)/(g.ssibry-g.ssimag);    
        tlim_interp = atan2(zinterp_lim-g.zmaxis,rinterp_lim-g.rmaxis);
        plot(psilim_interp,tlim_interp/pi,'bx','linewidth',3)
    
end





%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
% Read PSIN MIN FILE(S)
%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
max_file = ~isempty(fname_psimin2) + ~isempty(fname_psimin);
for ifile = 1:max_file
    if ifile == 1
        fname_psimin_tmp = fname_psimin;
    else
        fname_psimin_tmp = fname_psimin2;
    end    
    if ~isempty(fname_psimin_tmp)
        d2 = dlmread([run_path,fname_psimin_tmp]);
        psiNmin{ifile} = zeros(size(tline{ifile}));
        for i =1:numlines
            if d2(i+1) > 1000
               d2(i+1) = NaN;
            end
            psiNmin{ifile}(i,:) = d2(i+1);
            
        end
    end
    
    psiNmin_1d(ifile,:) = psiNmin{ifile}(:,1);
    psiN_1d(ifile,:) = psiNline{ifile}(:,1);
    r_1d(ifile,:) = rline{ifile}(:,1);
    z_1d(ifile,:) = zline{ifile}(:,1);
end


min1d = min(psiNmin_1d);
psi1d = psiN_1d(1,:);

figure; hold on; box on;
plot(psiN_1d(1,:),psiNmin_1d)
plot(psi1d,min1d,'k','linewidth',2)
plot([0,1],[0,1],'k--')
xlabel('\Psi_N','fontsize',14)
ylabel('\Psi_N^{min}','fontsize',14)

iend = find(isnan(min1d),1)-1;
iend2 = find(psi1d <= psiN_max_eval,1,'last');
iend = min([iend,iend2,length(psi1d)]);
m1 = min1d(1:iend);
p1 = psi1d(1:iend);


if 0  % -- should be obsolete compared to version below
    %------------------------------------------------------------------------------------------------------------------------------------
    %------------------------------------------------------------------------------------------------------------------------------------
    % Find island widths
    %------------------------------------------------------------------------------------------------------------------------------------
    %------------------------------------------------------------------------------------------------------------------------------------

    
    igo = 1;
    ind = 1;
    icount = 0;
    while igo == 1
        ifind = find(m1(ind+1:end) <= m1(ind),1);
        if ~isempty(ifind)
            icount = icount + 1;
            island_inds(icount,1:2) = [ind,ifind+ind];
            ind = ifind+ind;
        else
            ind = ind + 1;
        end
        if ind >= length(m1) - 1
            igo = 0;
        end
    end
    
    ni = size(island_inds,1);
    figure; hold on; box on;
    plot(p1,m1,'k.-','linewidth',2)
    for i = 1:ni
        plot(p1(island_inds(i,:)),m1(island_inds(i,:)),'r-','linewidth',2)
    end
    plot([0,1],[0,1],'k--')
    xlabel('\psi_N','fontsize',12)
    ylabel('\psi_N^{min}','fontsize',12)
    set(gca,'fontsize',12)
    title(mytitle)
end


%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
% Find island widths  ---> with interpolation
%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------



ninterp = 20000;
p2 = linspace(p1(1),p1(end),ninterp);
m2 = interp1(p1,m1,p2);

igo = 1;
ind = 1;
icount = 0;
while igo == 1
    ifind = find(m2(ind+1:end) <= m2(ind),1);
    if ~isempty(ifind)
        icount = icount + 1;
        island_inds2(icount,1:2) = [ind,ifind+ind];
        ind = ifind+ind + 1;
    else
        ind = ind + 1;
    end
    if ind >= length(m2) - 1
        igo = 0;
    end    
end

ni2 = size(island_inds2,1);
figure; hold on; box on;
% plot(p1,m1,'k-','linewidth',3)
plot(p2,m2,'k-','linewidth',2)
for i = 1:ni2
    plot(p2(island_inds2(i,:)),m2(island_inds2(i,:)),'r-','linewidth',2)
end
plot([0,1],[0,1],'k--')
psiN_islands = p2(island_inds2(:,:));
xlabel('\psi_N','fontsize',12)
ylabel('\psi_N^{min}','fontsize',12)
set(gca,'fontsize',12)
title(mytitle)

% figure; hold on; box on;
% plot(p2(1:end-1),diff(m2)./(p2(2)-p2(1)),'r')

%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
% Plot TS with island widths highlighted
%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
plotit = 0;
[R_ts,Z_ts,phi_ts,Te_ts,Te_std_ts,ne_ts,ne_std_ts,psiN_ts] = plot_ts_file_diiid(shot,times,plotit,g);
[R_ece,Z_ece,phi_ece,Te_ece,Te_std_ece,psiN_ece] = plot_ece_file_diiid(shot,times,plotit,g);
figure; hold on; box on;
subplot(2,1,1); hold on; box on;
if plot_ts
    errorbar(psiN_ts,Te_ts/1000,Te_std_ts/1000,'ko','linewidth',2)
end
if plot_ece
    errorbar(psiN_ece,Te_ece/1000,Te_std_ece/1000,'ko','linewidth',2)
end
title(['Shot: ',num2str(shot),', t=',num2str(times)],'fontsize',12)
ylabel('T_e (keV)','fontsize',12)
set(gca,'fontsize',12)
ylims = get(gca,'ylim');
for i=1:ni2
%     plot([1,1]*psiN_islands(i,1),ylims,'k')
%     plot([1,1]*psiN_islands(i,2),ylims,'k')
    patch([psiN_islands(i,1),psiN_islands(i,2),psiN_islands(i,2),psiN_islands(i,1),psiN_islands(i,1)],[ylims(1),ylims(1),ylims(2),ylims(2),ylims(1)],'b')
end
if plot_ts
    errorbar(psiN_ts,Te_ts/1000,Te_std_ts/1000,'ko','linewidth',2)
end
if plot_ece
    errorbar(psiN_ece,Te_ece/1000,Te_std_ece/1000,'ko','linewidth',2)
end

if plot_ts
    subplot(2,1,2); hold on; box on;
    errorbar(psiN_ts,ne_ts,ne_std_ts,'ko','linewidth',2)
    xlabel('\psi_N','fontsize',12);
    ylabel('n_e (10^{19} m^{-3})','fontsize',12)
    set(gca,'fontsize',12)
    ylims = get(gca,'ylim');
    for i=1:ni2
%         plot([1,1]*psiN_islands(i,1),ylims,'k')
%         plot([1,1]*psiN_islands(i,2),ylims,'k')
        patch([psiN_islands(i,1),psiN_islands(i,2),psiN_islands(i,2),psiN_islands(i,1),psiN_islands(i,1)],[ylims(1),ylims(1),ylims(2),ylims(2),ylims(1)],'b')
    end
    errorbar(psiN_ts,ne_ts,ne_std_ts,'ko','linewidth',2)
end



%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
% Plot poincare theta with islands
%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
figure; hold on; box on;
for ifile = 1:max_file
    plot(psiNline{ifile}.',tline{ifile}.'/pi,'k.','markersize',2)
end
xlabel('\psi_N','fontsize',12)
ylabel('\theta (\pi rad.)','fontsize',12)
set(gca,'fontsize',12)
title(mytitle)
% contour(g.r,g.z,psiN_g.',[0.5,0.6,0.7,0.8,0.9,1.0],'linewidth',2);

ylims = get(gca,'ylim');
for i=1:ni2
    plot([1,1]*psiN_islands(i,1),ylims,'c-','linewidth',2)
    plot([1,1]*psiN_islands(i,2),ylims,'c-','linewidth',2)
end

if plot_ts
    R(1:40)=1.94;
    Z=[0.0430,0.1720,0.2810,0.3320,0.3590,0.3910,0.4210,0.4770,0.5310,0.5630,0.5790,0.5950,0.6080,0.6230,0.6350,0.6410,0.6480,0.6540,0.6610,0.6680,...
        0.6740,0.6810,0.6870,0.6930,0.7010,0.7060,0.7130,0.7180,0.7250,0.7310,0.7370,0.7430,0.7490,0.7550,0.7630,0.7760,0.7880,0.8000,0.8130,0.8270];
    T = atan2(Z-g.zmaxis,R-g.rmaxis);
    PN = (g.ip_sign*get_psi_bicub(g,R,Z)-g.ssimag)/(g.ssibry-g.ssimag);
    plot(PN,T/pi,'bo-','linewidth',3)
    
    % Highlight along diagnostic line where islands are
    for i = 1:ni2
        pn_start = p2(island_inds2(i,1));
        pn_end = p2(island_inds2(i,2));
        Rh1 = interp1(PN,R,pn_start);
        Rh2 = interp1(PN,R,pn_end);
        Zh1 = interp1(PN,Z,pn_start);
        Zh2 = interp1(PN,Z,pn_end);        
        plot([pn_start,pn_end],[atan2(Zh1-g.zmaxis,Rh1-g.rmaxis),atan2(Zh2-g.zmaxis,Rh2-g.rmaxis)]./pi,'m-','linewidth',3)        
    end
    
end

if plot_ece
    Z2(1:40) = 0.004;
    R2 = [2.2399,2.2175,2.1951,2.1724,2.1496,2.1269,2.1045,2.0824,2.0608,2.0398,2.0193,1.9994,1.9801,1.9613,1.9430,1.9252,1.9252,1.9077,1.8905,1.8737,...
        1.8570,1.8406,1.8244,1.8084,1.7925,1.7768,1.7612,1.7458,1.7306,1.7154,1.7005,1.6856,1.6563,1.6276,1.5994,1.5717,1.5445,1.5178,1.4916,1.4658];
    T2 = atan2(Z2-g.zmaxis,R2-g.rmaxis);
    PN2 = (g.ip_sign*get_psi_bicub(g,R2,Z2)-g.ssimag)/(g.ssibry-g.ssimag);
%     [PN2,II2] = sort(PN2);
%     R2 = R2(II2);
%     Z2 = Z2(II2);
%     T2 = T2(II2);
    plot(PN2,T2/pi,'bo-','linewidth',3)

%     % Highlight along diagnostic line where islands are
%     for i = 1:ni2
%         pn_start = p2(island_inds2(i,1));
%         pn_end = p2(island_inds2(i,2));
%         Rh1 = interp1(PN2,R2,pn_start);
%         Rh2 = interp1(PN2,R2,pn_end);
%         Zh1 = interp1(PN2,Z2,pn_start);
%         Zh2 = interp1(PN2,Z2,pn_end);        
%         plot([pn_start,pn_end],[atan2(Zh1-g.zmaxis,Rh1-g.rmaxis),atan2(Zh2-g.zmaxis,Rh2-g.rmaxis)]./pi,'m-','linewidth',3)        
%     end
    
end


if plot_glim
    rlim = g.lim(1,:);
    zlim = g.lim(2,:);
    tlim = atan2(zlim-g.zmaxis,rlim-g.rmaxis);
    psilim = (g.ip_sign*get_psi_bicub(g,rlim,zlim)-g.ssimag)/(g.ssibry-g.ssimag);    
%     plot(psilim,tlim/pi,'c-','linewidth',3)
    
    ninterp_lim = 100;
    for i = 1:size(g.lim,2)-1
        rinterp_lim = linspace(g.lim(1,i),g.lim(1,i+1),ninterp_lim);
        zinterp_lim = linspace(g.lim(2,i),g.lim(2,i+1),ninterp_lim);
        psilim_interp = (g.ip_sign*get_psi_bicub(g,rinterp_lim,zinterp_lim)-g.ssimag)/(g.ssibry-g.ssimag);    
        tlim_interp = atan2(zinterp_lim-g.zmaxis,rinterp_lim-g.rmaxis);
        plot(psilim_interp,tlim_interp/pi,'bx','linewidth',3)
    end
        rinterp_lim = linspace(g.lim(1,end),g.lim(1,end),ninterp_lim);
        zinterp_lim = linspace(g.lim(2,1),g.lim(2,1),ninterp_lim);
        psilim_interp = (g.ip_sign*get_psi_bicub(g,rinterp_lim,zinterp_lim)-g.ssimag)/(g.ssibry-g.ssimag);    
        tlim_interp = atan2(zinterp_lim-g.zmaxis,rinterp_lim-g.rmaxis);
        plot(psilim_interp,tlim_interp/pi,'bx','linewidth',3)
end
    
for ifile = 1:max_file    
    clear poin_inds2
    poin_inds2 = round(interp1(psiNline{ifile}(:,1),1:length(psiNline{ifile}(:,1)),p2(island_inds2)));
    for i = 1:size(poin_inds2,1)*size(poin_inds2,2)
        plot(psiNline{ifile}(poin_inds2(i),:).',tline{ifile}(poin_inds2(i),:).'/pi,'r.')
    end
end
        




%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
% Plot poincare RZ with islands
%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
figure; hold on; box on;
for ifile = 1:max_file
    plot(rline{ifile},zline{ifile},'k.','markersize',2)
end
xlabel('R (m)','fontsize',12)
ylabel('Z (m)','fontsize',12)
set(gca,'fontsize',12)
% contour(g.r,g.z,psiN_g.',[0.5,0.6,0.7,0.8,0.9,1.0],'linewidth',2);
if plot_ts
    R(1:40)=1.94;
    Z=[0.0430,0.1720,0.2810,0.3320,0.3590,0.3910,0.4210,0.4770,0.5310,0.5630,0.5790,0.5950,0.6080,0.6230,0.6350,0.6410,0.6480,0.6540,0.6610,0.6680,...
        0.6740,0.6810,0.6870,0.6930,0.7010,0.7060,0.7130,0.7180,0.7250,0.7310,0.7370,0.7430,0.7490,0.7550,0.7630,0.7760,0.7880,0.8000,0.8130,0.8270];
    plot(R,Z,'ro','linewidth',2)
end
if plot_ece
    Z2(1:40) = 0.004;
    R2 = [2.2399,2.2175,2.1951,2.1724,2.1496,2.1269,2.1045,2.0824,2.0608,2.0398,2.0193,1.9994,1.9801,1.9613,1.9430,1.9252,1.9252,1.9077,1.8905,1.8737,...
        1.8570,1.8406,1.8244,1.8084,1.7925,1.7768,1.7612,1.7458,1.7306,1.7154,1.7005,1.6856,1.6563,1.6276,1.5994,1.5717,1.5445,1.5178,1.4916,1.4658];
    T2 = atan2(Z2-g.zmaxis,R2-g.rmaxis);
    PN2 = (g.ip_sign*get_psi_bicub(g,R2,Z2)-g.ssimag)/(g.ssibry-g.ssimag);
    plot(R2,Z2,'bo','linewidth',3)
end
if plot_glim
    plot(g.lim(1,:),g.lim(2,:),'k','linewidth',2)
end

for ifile = 1:max_file    
    clear poin_inds2
    poin_inds2 = round(interp1(psiNline{ifile}(:,1),1:length(psiNline{ifile}(:,1)),p2(island_inds2)));
    for i = 1:size(poin_inds2,1)*size(poin_inds2,2)
        plot(rline{ifile}(poin_inds2(i),:).',zline{ifile}(poin_inds2(i),:).','r.')
    end
end
        










% figure; hold on; box on;
% plot(psi1d,psi1d - min1d,'k','linewidth',2)

% figure; hold on; box on;
% for ifile = 1:max_file
%     plot(r_1d(ifile,:),z_1d(ifile,:),'ko')
% end

if 0

    
%     figure; hold on; box on;
%     contour(psiNline,tline/pi,psiNmin,100)

    p1 = linspace(0.4,1.0,100);
    t1 = linspace(-1,1.0,120);
    [p2,t2] = meshgrid(p1,t1);
    vq = griddata(psiNline(1:20000),tline(1:20000)/pi,psiNmin,p2,t2);
    contour(p2,t2,vq,100)

end

