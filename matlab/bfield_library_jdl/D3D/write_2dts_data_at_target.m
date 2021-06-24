clearvars;

% gfile_name = 'C:\Work\DIII-D\156855\g156855.04526_694';
gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Qprl experiment\power13mw\g174306.03500_152';
g = readg_g3d(gfile_name);

% switch 1
%     case 1
%         fname_out_part = '2MW_attach';
%         fname ='C:\Work\DIII-D\APS 2016\2dts_data\dts_fitted1d_osp_multishot_156855156856156857156858_remapto156855at4527msEFIT04_fuelD_ioncharge1_gamma7p539_see0p000_weightedfits.dat';        
%     case 2
%         fname_out_part = '2MW_detach';
%         fname ='C:\Work\DIII-D\APS 2016\2dts_data\dts_fitted1d_osp_multishot_156859156861156862156863_remapto156855at4527msEFIT04_fuelD_ioncharge1_gamma7p539_see0p000_weightedfits.dat';
%     case 3
%         fname_out_part = '5MW_attach';
%         fname ='C:\Work\DIII-D\APS 2016\2dts_data\dts_fitted1d_osp_multishot_156866156867156868156869_remapto156855at4527msEFIT04_fuelD_ioncharge1_gamma7p539_see0p000_weightedfits.dat';
%     case 4
%         fname_out_part = '5MW_detach';
%         fname ='C:\Work\DIII-D\APS 2016\2dts_data\dts_fitted1d_osp_multishot_156871156872156873156874_remapto156855at4527msEFIT04_fuelD_ioncharge1_gamma7p539_see0p000_weightedfits.dat';        
% end
run_path = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Qprl experiment\1743XX\';

fname = fullfile(run_path,'dts_fitted1d_osp_174306_2000to5000msEFIT03_remapto3500msEFIT03_fuelD_ioncharge1_gamma7p539_see0p000_weightedfits.dat');
outfile = fullfile(run_path,'dts_osp_174306.mat');

% fname = fullfile(run_path,'dts_fitted1d_osp_174310_3000to4800msEFIT03_remapto174306at3500msEFIT03_fuelD_ioncharge1_gamma7p539_see0p000_weightedfits.dat');
% outfile = fullfile(run_path,'dts_osp_174310.mat');


ts2d = read_2dts_dat_file(fname);

WRITE_IT = 0;

plot_log = 0;

minval_Te = 0.5;
maxval_Te = 200;

minval_ne = 0.02;  % e 20
maxval_ne = 1;

minval_q = 0;
maxval_q = 10;

PLOT_TE =1;
PLOT_NE = 1;
PLOT_Q = 1;

% colorMapName = 'parula';
colorMapName = 'plasma';
% colorMapName = 'rainbow';

% ichan_use = 0;  %0:7
ichan_use = 0:7;  %0:7

% Make 1D arrays
icount = 1;
for ichan = ichan_use
    for i = 1:length(ts2d.chan{ichan+1}.Te)
        if ts2d.chan{ichan+1}.Rremap(i) > 1e-3
            R1d(icount) = ts2d.chan{ichan+1}.Rremap(i);
            R1d_maptarg(icount) = ts2d.chan{ichan+1}.Rmap_targ(i);
            Z1d(icount) = ts2d.chan{ichan+1}.Zremap(i);
            Z1d_maptarg(icount) = ts2d.chan{ichan+1}.Zmap_targ(i);
            Te1d(icount) = ts2d.chan{ichan+1}.Te(i);
            dTe1d(icount) = ts2d.chan{ichan+1}.dTe(i);
            ne1d(icount) = ts2d.chan{ichan+1}.ne(i);
            dne1d(icount) = ts2d.chan{ichan+1}.dne(i);
            q1d(icount) = ts2d.chan{ichan+1}.qperp(i);
            dq1d(icount) = ts2d.chan{ichan+1}.qperp_err(i);
            time1d(icount) = ts2d.chan{ichan+1}.time(i);
            icount = icount + 1;
        end
    end
end

if PLOT_TE
    figure; hold on; box on;
    if plot_log
        scatter(R1d,Z1d,[],log10(Te1d),'filled')
    else
        scatter(R1d,Z1d,[],Te1d,'filled')
    end
    colorbar;
    ncols = 1024;
    colormap(colorflipper(ncols,colorMapName));
    if plot_log
        set(gca,'clim',log10([minval_Te,maxval_Te]))
    else
        set(gca,'clim',[minval_Te,maxval_Te])
    end
    plot_sep_g(g,1,0);
    axis([1.35,1.65,-1.25,-0.95])
    title('T_e [eV]')
    
    

end

if PLOT_NE
    figure; hold on; box on;
    if plot_log
        scatter(R1d,Z1d,[],log10(ne1d),'filled')
    else
        scatter(R1d,Z1d,[],ne1d,'filled')
    end
    colorbar;
    ncols = 1024;
    colormap(colorflipper(ncols,colorMapName));
    if plot_log
        set(gca,'clim',log10([minval_ne,maxval_ne]))
    else
        set(gca,'clim',[minval_ne,maxval_ne])
    end
    plot_sep_g(g,1,0);
    axis([1.35,1.65,-1.25,-0.95])
    title('n_e [10^{20} m^{-3}]')
end

if PLOT_Q
    figure; hold on; box on;
    if plot_log
        scatter(R1d,Z1d,[],log10(q1d),'filled')
    else
        scatter(R1d,Z1d,[],q1d,'filled')
    end
    colorbar;
    ncols = 1024;
    colormap(colorflipper(ncols,colorMapName));
    if plot_log
        set(gca,'clim',log10([minval_q,maxval_q]))
    else
        set(gca,'clim',[minval_q,maxval_q])
    end    
    plot_sep_g(g,1,0);
    axis([1.35,1.65,-1.25,-0.95])
    title('q perp [MW/m^2]')
end


if length(ichan_use) == 1
    if ichan_use == 0
        
        %         figure; hold on; box on;
        %         errorbar(R1d_maptarg,Te1d,dTe1d,'rx')
        %         errorbar(R1d,Te1d,dTe1d,'ko')
        %         legend('Rmap targ','Rremap')
        %         xlabel('R [m]')
        %         ylabel('T_e [eV]')
        %
        %         figure; hold on; box on;
        %         errorbar(R1d_maptarg,ne1d,dne1d,'rx')
        %         errorbar(R1d,ne1d,dne1d,'ko')
        %         legend('Rmap targ','Rremap')
        %         xlabel('R [m]')
        %         ylabel('n_e [10^{20} m^{-3}]')
        %
        %         figure; hold on; box on;
        %         errorbar(R1d_maptarg,q1d,dq1d,'rx')
        %         errorbar(R1d,q1d,dq1d,'ko')
        %         legend('Rmap targ','Rremap')
        %         xlabel('R [m]')
        %         ylabel('q_{\perp} [MW/m^2]')
        
        
        
        
        Rsep = interp1(ts2d.chan{1}.psiN,ts2d.chan{1}.Rmap_targ,1);
        figure; hold on; box on;
        errorbar((R1d-Rsep)*100,q1d,dq1d,'ko')
        xlabel('R [m]')
        ylabel('q_{\perp} [MW/m^2]')
        
        figure; hold on; box on;
        errorbar((R1d-Rsep)*100,Te1d,dTe1d,'ko')
        xlabel('R [m]')
        ylabel('T_{e} [eV]')
        
        figure; hold on; box on;
        errorbar((R1d-Rsep)*100,ne1d,dne1d,'ko')
        xlabel('R [m]')
        ylabel('n_{e} []')        
        
        % experimental
        figure; hold on; box on;
        scatter((R1d-Rsep)*100,q1d,[],time1d,'filled')
        colorbar;
        ncols = 1024;
        colormap(colorflipper(ncols,colorMapName));
        
      
    end
end



if WRITE_IT
    [pathstr] = fileparts(fname);
    
    save(outfile,'ts2d');
    
    
    if length(ichan_use) == 1 & 0 % Switched to mat version
        
        if ichan_use == 0
                        
            fname_out = fullfile(pathstr,strcat('Te_2DTS_1D_chan0_targ_',fname_out_part,'.dat'));
            fprintf('Writing Te data to file: %s\n',fname_out)
            fid = fopen(fname_out,'w');
            fprintf(fid,'%d\n',length(Te1d));
            for i = 1:length(Te1d)
                fprintf(fid,'%e %e %e\n',R1d_maptarg(i),Te1d(i),dTe1d(i));
            end
            fclose(fid);
            
            fname_out = fullfile(pathstr,strcat('ne_2DTS_1D_chan0_targ_',fname_out_part,'.dat'));
            fprintf('Writing ne data to file: %s\n',fname_out)
            fid = fopen(fname_out,'w');
            fprintf(fid,'%d\n',length(ne1d));
            for i = 1:length(ne1d)
                fprintf(fid,'%e %e %e\n',R1d_maptarg(i),ne1d(i),dne1d(i));
            end
            fclose(fid);
            
            fname_out = fullfile(pathstr,strcat('qperp_2DTS_1D_chan0_targ_',fname_out_part,'.dat'));
            fprintf('Writing ne data to file: %s\n',fname_out)
            fid = fopen(fname_out,'w');
            fprintf(fid,'%d\n',length(q1d));
            for i = 1:length(q1d)
                fprintf(fid,'%e %e %e\n',R1d_maptarg(i),q1d(i),dq1d(i));
            end
            fclose(fid);
                     
        else
            fprintf('Skipping write because not channel 0!\n')
            
            
        end
    else
        fprintf('Skipping write because too many channels!\n')
    end
end


