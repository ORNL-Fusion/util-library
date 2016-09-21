function sep = plot_sep_g(g,plotit,newfig)
if nargin < 2
    plotit = 1;
end
if nargin < 3
    newfig = 0;
end

%(run_path,g,newfig,plotit,sep_file_name,rguess,zguess,de,SECOND_SEP,save_in_dir)



clip_at_glim = 1;

if newfig == 1
    figure; hold on;   
end



box_r = [g.r(1),g.r(g.mw),g.r(g.mw),g.r(1),g.r(1)];
box_z = [g.z(1),g.z(1),g.z(g.mh),g.z(g.mh),g.z(1)];

[pathstr,fn,ext] = fileparts(g.filename);
fname = fullfile(pathstr,strcat(fn,ext,'_sep'));

fid = fopen(fname,'r');
if fid == -1
    disp(['Did not find sep. data file: ',fname])
    disp('Generating separatrix file')
    
    lim_r = g.lim(1,g.lim(1,:) > 0);
    lim_z = g.lim(2,g.lim(1,:) > 0);
    
%     plot(box_r,box_z,'c-')
%     plot(lim_r,lim_z,'r-')
    xpt_info = find_xpt_jl(g,1,1,1e-8,0);
    xr1 = xpt_info.rx;
    xz1 = xpt_info.zx;
    
    
    dxr = 1e-3;
    dphi_fl = g.ip_sign*0.1*pi/180;
    ns_fl = floor(0.3*pi/abs(dphi_fl));  % Factor of 0.x to stop from leaving grid
    maxtries = 1000;
    
    for idir = [1,-1]
        disp(['Working on direction: ',num2str(idir)])
        
        rstart = xr1-dxr;
        zstart = xz1;
        phistart = 0;

        rsep0 = [];
        zsep0 = [];
        for ii = 1:maxtries
            bfield.type = 'gfile';
            bfield.g = g;
            [stmp0,ierr_ff]=follow_fieldlines_rzphi_dphi(bfield,rstart,zstart,phistart,idir*dphi_fl,ns_fl);
    
            rsep0 = [rsep0;stmp0.r(2:end)];
            zsep0 = [zsep0;stmp0.z(2:end)];
            rstart = stmp0.r(end);
            zstart = stmp0.z(end);

            inp1 = inpolygon(rstart,zstart,lim_r,lim_z);
            inp2 = inpolygon(rstart,zstart,box_r,box_z);
            if inp1 == 0 || inp2 == 0
                break
            end
        end
        if ii == maxtries
            error('increase maxtries')
        end        
%         plot(rsep0,zsep0,'g')
        if idir == 1
            rsep1 = rsep0;
            zsep1 = zsep0;
        end
    end
    rsep2=[flipud(rsep0(2:end));rsep1];
    zsep2=[flipud(zsep0(2:end));zsep1];

    nsep_want = 1000;
    rsep = rsep2(1:floor(length(rsep2)/nsep_want):end);
    zsep = zsep2(1:floor(length(zsep2)/nsep_want):end);
    
    if clip_at_glim
        inp = inpolygon(rsep,zsep,g.lim(1,:),g.lim(2,:));
        rsep = rsep(inp); zsep = zsep(inp); nsep = length(rsep);
    end

    if plotit == 1
        plot(rsep,zsep,'k-')
    end
    disp('>>>> Writing sep file')
    fid = fopen(fname,'w');
    fprintf(fid,'%i\n',length(rsep));
    fprintf(fid,'%f ',rsep);
    fprintf(fid,'\n');
    fprintf(fid,'%f ',zsep);
    fprintf(fid,'\n');
    fclose(fid);
    nsep = length(rsep);
    
    sep.rsep = rsep;
    sep.zsep = zsep;
    sep.nsep = nsep;
    
    
else
    nsep = fscanf(fid,'%i\n',1);
    rsep = fscanf(fid,'%f\n',nsep);
    zsep = fscanf(fid,'%f\n',nsep);
    fclose(fid);
    if plotit ==1
        plot(rsep,zsep,'k-')
    end
    sep.rsep = rsep;
    sep.zsep = zsep;
    sep.nsep = nsep;
        
end
if nargout == 0
    clear sep
end
