function ts2d = read_2dts_dat_file(fname)
% CHANNEL INDICES ARE OFFSET BY 1!!!!
% fname = 'C:\Work\DIII-D\164723\2DTS_Stuff\dts_fitted1d_osp_164722_2600to5000msEFIT01_remapto141905at3500msEFIT01_fuelD_ioncharge1_gamma7p539_see0p000.dat';

% fname ='C:\Work\DIII-D\APS 2016\2dts_data\dts_fitted1d_osp_multishot_156855156856156857156858_remapto156855at4527msEFIT04_fuelD_ioncharge1_gamma7p539_see0p000_weightedfits.dat';
fid = fopen(fname);

fgetl(fid); % Title
fgetl(fid); % blank
labels = fgetl(fid); % legend


icount = 0;
iline = 3;

wid = 20; % Width of each entry
% ncol = 86;


while ~feof(fid)
    
    line = fgetl(fid); iline = iline + 1;
    if isempty(line)
        line = fgetl(fid); iline = iline + 1;  % Probably blank line between channels
    end
    ncol = length(line)/wid;
    if ncol == 86
        ver = 1;
        iwants = [1:3,8,24:29,30:32,35:37,40,80:81];  % Select columns to be used, must be recognized below
    elseif ncol == 84
        ver = 2;
        iwants = [1:3,8,24:29,30:32,35:37,40];  % Select columns to be used, must be recognized below
    elseif ncol == 91
        ver = 3;
        iwants = [1:3,8,24:29,30:32,35:37,40];  % Select columns to be used, must be recognized below
    else
         if feof(fid)
            break;
         end
         fprintf('ncol found %d\n',ncol)
        error('Do not understand file');
    end

    
    icount = icount + 1;
    for i = iwants
        inds = (i-1)*wid+1:i*wid;
        mylabel = labels(inds);
        if iline == 4
            fprintf('Returning data from column with legend: %s\n',mylabel)
        end
        switch i
            case 1  % Shot
                shot(icount) = sscanf(line(inds),'%d',1);
            case 2  % Laser time
                time(icount) = sscanf(line(inds),'%f',1);
            case 3  % Channel(s)
                num_channels = length(find(line(inds) == ',')) + 1;
                myfmt = '%d';
                for j = 2:num_channels
                    myfmt = strcat(myfmt,',%d');
                end
                channels{icount} = sscanf(line(inds),myfmt,num_channels);
            case 8   % Psi_N
                PsiN(icount) = sscanf(line(inds),'%f',1);
            case 24  % R_remap [m]
                R_remap_m(icount) = sscanf(line(inds),'%f',1);
            case 25  % Z_remap [m]
                Z_remap_m(icount) = sscanf(line(inds),'%f',1);
            case 26  % R_map,targ [m]
                R_map_targ_m(icount) = sscanf(line(inds),'%f',1);
            case 27  % Z_map,targ [m]
                Z_map_targ_m(icount) = sscanf(line(inds),'%f',1);                
            case 28  % R-Rsep,map [m]
                RmRsep_map_targ_m(icount) = sscanf(line(inds),'%f',1);
            case 29  % Z-Zsep,map [m]
                ZmZsep_map_targ_m(icount) = sscanf(line(inds),'%f',1);                                
            case 30  % R_OMP [m]
                Romp_m(icount) = sscanf(line(inds),'%f',1);
            case 31  % Te,meas [eV]
                Te_eV(icount) = sscanf(line(inds),'%f',1);
            case 32  % Te,meas,err [eV]
                Te_err_eV(icount) = sscanf(line(inds),'%f',1);
            case 35  % Te,fit [eV]
                Te_fit_eV(icount) = sscanf(line(inds),'%f',1);
            case 36  % ne,meas [e20,m^-3]
                ne_e20(icount) = sscanf(line(inds),'%f',1);
            case 37  % ne,err [e20,m^-3]
                ne_err_e20(icount) = sscanf(line(inds),'%f',1);
            case 40  % ne,fit [e20,m^-3]
                ne_fit_e20(icount) = sscanf(line(inds),'%f',1);
            case 80  % qperp (MW/m^2)  
                qperp_MWm2(icount) = sscanf(line(inds),'%f',1);
            case 81  % qperp,err (MW/m^2)
                qperp_err_MWm2(icount) = sscanf(line(inds),'%f',1);                
            otherwise
                warning('No case set up for i=%d\n',i)
        end
    end         
    
end

fclose(fid);


% Just use remaps to one channel
icount_chans = zeros(8,1);
i2 = 1;
for i = 1:length(channels)
    if length(channels{i}) == 1
        ic = channels{i}(1)+1;
        icount_chans(ic) = icount_chans(ic) + 1;
        ts2d.channel_array(i2)                    = channels{i}(1);
        ts2d.chan{ic}.shot(icount_chans(ic))      = shot(i);
        ts2d.chan{ic}.time(icount_chans(ic))      = time(i);
        ts2d.chan{ic}.Te(icount_chans(ic))        = Te_eV(i);
        ts2d.chan{ic}.dTe(icount_chans(ic))       = Te_err_eV(i);
        ts2d.chan{ic}.Tefit(icount_chans(ic))     =  Te_fit_eV(i);
        ts2d.chan{ic}.ne(icount_chans(ic))        = ne_e20(i);
        ts2d.chan{ic}.dne(icount_chans(ic))       = ne_err_e20(i);
        ts2d.chan{ic}.nefit(icount_chans(ic))     = ne_fit_e20(i);
        ts2d.chan{ic}.Romp(icount_chans(ic))      = Romp_m(i);
        ts2d.chan{ic}.psiN(icount_chans(ic))      = PsiN(i);
        ts2d.chan{ic}.Rremap(icount_chans(ic))    = R_remap_m(i);
        ts2d.chan{ic}.Zremap(icount_chans(ic))    = Z_remap_m(i);
        ts2d.chan{ic}.Rmap_targ(icount_chans(ic)) = R_map_targ_m(i);
        ts2d.chan{ic}.Zmap_targ(icount_chans(ic)) = Z_map_targ_m(i);
        ts2d.chan{ic}.RmRsep_map(icount_chans(ic)) = RmRsep_map_targ_m(i);
        ts2d.chan{ic}.ZmZsep_map(icount_chans(ic)) = ZmZsep_map_targ_m(i);
        
        if ver == 1
        ts2d.chan{ic}.qperp(icount_chans(ic))     = qperp_MWm2(i);
        ts2d.chan{ic}.qperp_err(icount_chans(ic)) = qperp_err_MWm2(i);
        end
        i2 = i2 + 1;
    end
end


