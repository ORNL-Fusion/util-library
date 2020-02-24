clearvars;

plot_ts = 0;
plot_ece = 0;
plot_glim = 0;

fname_psimin = []; fname_psimin2 = []; gfile_name =[];

if 1
%    run_path = '/home/jjl/util-library/fortran/poincare/examples/M3DC1_14170'; suffix = 'test'; mytitle = 'test';
   run_path = '/home/jjl/util-library/fortran/poincare/examples/M3DC1_14170'; suffix = 'm3dc1_test'; mytitle = 'test';
end

if 0
    % 164723 <------------------------*********************--------------------
    gfile_name = 'C:\Work\DIII-D\164723\g164723.03059_410';
    % run_path = 'C:\Work\fortran\test_poincare\164723\1FL_AS\'; suffix = '30deg2'; mytitle = '1FL \phi = 30';
    %  run_path = 'C:\Work\fortran\test_poincare\164723\1FL\'; suffix = '30deg'; mytitle = '1FL \phi = 30';
    % run_path = 'C:\Work\fortran\test_poincare\164723\2FL\'; suffix = '30deg'; mytitle = '2FL \phi = 30';
    % run_path = 'C:\Work\fortran\test_poincare\164723\AS\'; suffix = '0deg'; mytitle = 'AS \phi = 0';
    run_path = 'C:\Work\fortran\test_poincare\164723\VAC\'; suffix = '0deg';  mytitle = 'VAC \phi = 0';
    % run_path = 'C:\Work\fortran\test_poincare\164723\1FL_VAC\'; suffix = '30deg';  mytitle = '1FL VAC \phi = 30';
    % run_path = 'C:\Work\fortran\test_poincare\164723\ipec_high_eq\'; suffix = '0deg'; mytitle = 'IPEC eq \phi = 0';
    % run_path = 'C:\Work\fortran\test_poincare\164723\ipec_high_vac\'; suffix = '30deg'; mytitle = 'IPEC VAC \phi = 30';
    % run_path = 'C:\Work\fortran\test_poincare\164723\ipec_high_pert\'; suffix = '0deg'; mytitle = 'IPEC pert \phi = 0';
    % run_path = 'C:\Work\fortran\test_poincare\164723\gpec_high_pert\'; suffix = '0deg'; mytitle = 'GPEC pert \phi = 0';
    % run_path = 'C:\Work\fortran\test_poincare\164723\XPAND_pert\'; suffix = '0deg'; mytitle = 'XPAND pert \phi = 0';
    % run_path = 'C:\Work\fortran\test_poincare\164723\XPAND_vac\'; suffix = '0deg'; mytitle = 'XPAND vac \phi = 0';
end

if 0
    % 160884 <------------------------*********************--------------------
    % gfile_name = 'C:\Work\DIII-D\160884\efits\g160884.05009_537';
    % run_path = 'C:\Work\fortran\test_poincare\160884\5009\VAC\'; suffix = '0deg';  mytitle = 'VAC \phi = 0';
    gfile_name = 'C:\Work\M3DC1\160884\5000\g160884.05000_m3dc1';
    % run_path = 'C:\Work\fortran\test_poincare\160884\5009\1FL\'; suffix = '30deg';  mytitle = '1FL \phi = 30';
    % run_path = 'C:\Work\fortran\test_poincare\160884\5009\VAC_m3dc1\'; suffix = '0deg';  mytitle = 'VAC_m3dc1 \phi = 0';
    % run_path = 'C:\Work\fortran\test_poincare\160884\5009\VAC_m3dc1\'; suffix = '30deg';  mytitle = 'VAC_m3dc1 \phi = 30';
    run_path = 'C:\Work\fortran\test_poincare\160884\5009\VAC_m3dc1\'; suffix = '-30deg';  mytitle = 'VAC_m3dc1 \phi = -30';
end

% run_path = 'C:\Work\fortran\test_poincare\165274\'; suffix = '-120'; gfile_name = 'C:\Work\DIII-D\165274\kinetic\g165274.02120'; mytitle = '\phi = 0'; plot_ts = 0; plot_ece=0;

if 0
    % 154929 <------------------------*********************--------------------
    run_path = 'C:\Work\fortran\test_poincare\154929\rmp\'; suffix = '_60_n100_fake_diag2'; mytitle='test';shot = 154921; times = 3750;
    %     run_path = 'C:\Work\fortran\test_poincare\154929\rmp\'; suffix = '240_n250_fake_diag2'; mytitle='test';shot = 154921; times = 3750;
    gfile_name = 'C:\Work\DIII-D\154929\kinetic\g154929.03750';
end

if 0
    run_path = 'C:\Work\fortran\test_poincare\vmec_coils_to_fil\'; suffix = ''; mytitle = 'test';
end

if 0
    % 165274 <------------------------*********************--------------------
    %     run_path = 'C:\Work\fortran\test_poincare\165274\'; suffix = 'test'; mytitle='165274';shot = 165274;% times = 3750;
    run_path = 'C:\Work\fortran\test_poincare\165274\'; suffix = '-120'; mytitle='165274';shot = 165274;% times = 3750;
    %     run_path = 'C:\Work\fortran\test_poincare\165274\'; suffix = '2_1'; mytitle='2/1';shot = 165274;% times = 3750;
    %     run_path = 'C:\Work\fortran\test_poincare\165274\'; suffix = '3_1'; mytitle='3/1';shot = 165274;% times = 3750;
    %     run_path = 'C:\Work\fortran\test_poincare\165274\'; suffix = '4_1'; mytitle='4/1';shot = 165274;% times = 3750;
    %     run_path = 'C:\Work\fortran\test_poincare\165274\'; suffix = '4_1_3turns'; mytitle='4/1';shot = 165274;% times = 3750;
    %     run_path = 'C:\Work\fortran\test_poincare\165274\'; suffix = 'normal_to_compare_to_3'; mytitle='4/1';shot = 165274;% times = 3750;
    %     run_path = 'C:\Work\fortran\test_poincare\165274\'; suffix = '3turn_test'; mytitle='4/1';shot = 165274;% times = 3750;
    gfile_name = 'C:\Work\DIII-D\165274\kinetic\g165274.02120';
end

if 0
    %     run_path = 'C:\Work\fortran\test_poincare\W7-X\OP2\22kA\bmw\'; suffix = '60surf';mytitle = 'bmw';
    % run_path = 'C:\Work\fortran\test_poincare\W7-X\OP2\22kA\mfbe\'; suffix = '60surf';mytitle = 'mfbe';
    % run_path = 'C:\Work\fortran\test_poincare\W7-X\OP2\22kA\extender\'; suffix = '40surf_v2';mytitle = 'extender';
    % run_path = 'C:\Work_archive\fortran\test_poincare\W7-X\OP2\22kA\extender\'; suffix = '60surf';mytitle = 'extender';
    % run_path = 'C:\Work\fortran\test_poincare\W7-X\OP2\22kA\bmw\'; suffix = '40surf_18deg';mytitle = 'bmw';
    % run_path = 'C:\Work_archive\fortran\test_poincare\W7-X\OP2\22kA_altern\'; suffix = '0deg_200tran_60surf';mytitle = 'extender';
    % run_path = 'C:\Work_archive\fortran\test_poincare\W7-X\OP12a\22kA_mimic\'; suffix = '0deg_200tran_60surf';mytitle = 'extender';% 1.75 parts
    % run_path = 'C:\Work_archive\fortran\test_poincare\W7-X\OP12a\22kA_mimic\'; suffix = '-160.75_200tran_60surf';mytitle = 'extender'; % 16.75 parts
    
    % MPM scan
    run_path = 'C:\Work_archive\fortran\test_poincare\W7-X\OP12a\0kA_mimic\'; suffix = '15.2deg_test1';mytitle = '';
    % run_path = 'C:\Work_archive\fortran\test_poincare\W7-X\OP12a\11kA_mimic\'; suffix = '15.2deg_test1';mytitle = '11kA';
    % run_path = 'C:\Work_archive\fortran\test_poincare\W7-X\OP12a\22kA_mimic\'; suffix = '15.2deg_test1';mytitle = '';
    % run_path = 'C:\Work_archive\fortran\test_poincare\W7-X\OP12a\32kA_mimic\'; suffix = '15.2deg_test1';mytitle = '32kA';
    % run_path = 'C:\Work_archive\fortran\test_poincare\W7-X\OP12a\43kA_mimic\'; suffix = '15.2deg_test1';mytitle = '43kA';
    
    % run_path = 'C:\Work_archive\fortran\test_poincare\W7-X\OP12b\43kA_Icc_scan\Icc_+2500\'; suffix = '0deg_200tran_100surf'; mytitle = '0deg +2500';
    % run_path = 'C:\Work_archive\fortran\test_poincare\W7-X\OP12b\43kA_Icc_scan\Icc_+1250\'; suffix = '0deg_200tran_100surf'; mytitle = '0deg +1250';
    % run_path = 'C:\Work_archive\fortran\test_poincare\W7-X\OP12b\43kA_Icc_scan\Icc_0\'; suffix = '0deg_200tran_100surf'; mytitle = '0deg 0';
    % run_path = 'C:\Work_archive\fortran\test_poincare\W7-X\OP12b\43kA_Icc_scan\Icc_-1250\'; suffix = '0deg_200tran_100surf'; mytitle = '0deg -1250';
    % run_path = 'C:\Work_archive\fortran\test_poincare\W7-X\OP12b\43kA_Icc_scan\Icc_-2500\'; suffix = '0deg_200tran_100surf'; mytitle = '0deg -2500';
    
    % run_path = 'C:\Work_archive\fortran\test_poincare\W7-X\OP2\22kA_altern\'; suffix = '-145.75deg_200tran_60surf';mytitle = 'extender';
end

if 0
    %     run_path = 'C:\Work\fortran\test_poincare\W7-X\OP2\22kA\bmw\'; suffix = 'find_lcfs';mytitle = 'bmw';
    %     run_path = 'C:\Work\fortran\test_poincare\W7-X\OP2\22kA\mfbe\'; suffix = 'find_lcfs';mytitle = 'mfbe';
end

MYCOL = 'k';
NEWFIG = 1;
TITLE = 1;
AXES_LABELS = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% DONE SETTING RUN PATH
%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(suffix)
    fname_psimin = ['psiN_min_output_',suffix,'.out'];fname_psimin2 = ['psiN_min_output2_',suffix,'.out'];
    fname = ['poincare_output_',suffix,'.out'];  fname2 = ['poincare_output2_',suffix,'.out'];
else
    fname_psimin = 'psiN_min_output.out';fname_psimin2 = 'psiN_min_output2.out';
    fname = 'poincare_output.out';  fname2 = 'poincare_output2.out';
end

psiN_max_eval = 0.98;

if ~isempty(gfile_name)
    g=readg_g3d(gfile_name);
    psiN_g = (g.psirz-g.ssimag)/(g.ssibry-g.ssimag);
else
    g = [];
end


poinc = read_poincare_file(run_path,fname,fname2,g);

%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
% PLOT RZ poincare
%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------

if 1
    if NEWFIG; figure; hold on; box on; end
    for ifile = 1:2
        plot(poinc.rline{ifile},poinc.zline{ifile},strcat(MYCOL,'.'),'markersize',2)
    end
    if AXES_LABELS
        xlabel('R (m)','fontsize',12)
        ylabel('Z (m)','fontsize',12)
        set(gca,'fontsize',12)
    end
    if TITLE
        title(mytitle)
    end
    %     contour(g.r,g.z,psiN_g.',[0.5,0.6,0.7,0.8,0.9,1.0],'linewidth',2);
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
if 0
    
    for ifile = 1:2
        tline{ifile} = atan2(poinc.zline{ifile}-g.zmaxis,poinc.rline{ifile}-g.rmaxis);
    end
    
    
    figure; hold on; box on;
    for ifile = 1:2
        plot(poinc.psiNline{ifile}.',tline{ifile}.'/pi,'k.','markersize',2)
    end
    xlabel('\psi_N','fontsize',12)
    ylabel('\theta (\pi rad.)','fontsize',12)
    set(gca,'fontsize',12)
    title(mytitle)
    % contour(g.r,g.z,psiN_g.',[0.5,0.6,0.7,0.8,0.9,1.0],'linewidth',2);
    
    
    %     figure; hold on; box on;
    %     c = colormap(colorflipper(256,'jet'));
    %     for ifile = 1:2
    %         for i = 1:size(tline{ifile},1)
    %             cind = round(interp1([0,1],[1,size(c,1)],i/(size(tline{ifile},1)+1)));
    %             plot(poinc.psiNline{ifile}(i,:),tline{ifile}(i,:)/pi,'.','markersize',2,'color',c(cind,:))
    %         end
    %     end
    %     xlabel('\psi_N','fontsize',12)
    %     ylabel('\theta (\pi rad.)','fontsize',12)
    %     set(gca,'fontsize',12)
    %     title(mytitle)
    
    
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
    
    
end


%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
% Read PSIN MIN FILE(S)
%------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------
if 0
    %     2 = ~isempty(fname_psimin2) + ~isempty(fname_psimin);
    for ifile = 1:2
        if ifile == 1
            fname_psimin_tmp = fname_psimin;
        else
            fname_psimin_tmp = fname_psimin2;
        end
        if ~isempty(fname_psimin_tmp)
            d2 = dlmread([run_path,fname_psimin_tmp]);
            psiNmin{ifile} = zeros(size(tline{ifile}));
            for i =1:poinc.numlines
                if d2(i+1) > 1000
                    d2(i+1) = NaN;
                end
                psiNmin{ifile}(i,:) = d2(i+1);
                
            end
        end
        
        psiNmin_1d(ifile,:) = psiNmin{ifile}(:,1);
        psiN_1d(ifile,:) = poinc.psiNline{ifile}(:,1);
        r_1d(ifile,:) = poinc.rline{ifile}(:,1);
        z_1d(ifile,:) = poinc.zline{ifile}(:,1);
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
        errorbar(psiN_ts,Te_ts,Te_std_ts,'ko','linewidth',2)
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
    for ifile = 1:2
        plot(poinc.psiNline{ifile}.',tline{ifile}.'/pi,'k.','markersize',2)
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
    
    for ifile = 1:2
        clear poin_inds2
        poin_inds2 = round(interp1(poinc.psiNline{ifile}(:,1),1:length(poinc.psiNline{ifile}(:,1)),p2(island_inds2)));
        for i = 1:size(poin_inds2,1)*size(poin_inds2,2)
            plot(poinc.psiNline{ifile}(poin_inds2(i),:).',tline{ifile}(poin_inds2(i),:).'/pi,'r.')
        end
    end
    
    
    
    
    
    %------------------------------------------------------------------------------------------------------------------------------------
    %------------------------------------------------------------------------------------------------------------------------------------
    % Plot poincare RZ with islands
    %------------------------------------------------------------------------------------------------------------------------------------
    %------------------------------------------------------------------------------------------------------------------------------------
    figure; hold on; box on;
    for ifile = 1:2
        plot(poinc.rline{ifile},poinc.zline{ifile},'k.','markersize',2)
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
    
    for ifile = 1:2
        clear poin_inds2
        poin_inds2 = round(interp1(poinc.psiNline{ifile}(:,1),1:length(poinc.psiNline{ifile}(:,1)),p2(island_inds2)));
        for i = 1:size(poin_inds2,1)*size(poin_inds2,2)
            plot(poinc.rline{ifile}(poin_inds2(i),:).',poinc.zline{ifile}(poin_inds2(i),:).','r.')
        end
    end
    
    
    
    
    
    
end




% figure; hold on; box on;
% plot(psi1d,psi1d - min1d,'k','linewidth',2)

% figure; hold on; box on;
% for ifile = 1:2
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

