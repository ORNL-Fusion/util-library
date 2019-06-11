clearvars;

ir_path = 'C:\Users\jjl\Dropbox (ORNL)\W7X\OP12a_Data\Jeremy_H';
scene_path = ir_path;

PLOT_UNAVERAGED = 0;

% myshots = [41,39,32,44,21];  % H 2MW density scan (0,11,22,32,43)
% myshots = [14,15,16];
myshots = [32];

CAM = 4;

if CAM == 4
    fname_scene_h = fullfile(scene_path,'s4lh_3Dlines_coordinate.txt');
    fname_scene_v = fullfile(scene_path,'s4lv_3Dlines_coordinate.txt');
elseif CAM == 2
    fname_scene_h = fullfile(scene_path,'s2uh_3D.txt');
    fname_scene_v = fullfile(scene_path,'s2uv_3D.txt');
elseif CAM == 3
    fname_scene_h = fullfile(scene_path,'s3lh_3D.txt');
    fname_scene_v = fullfile(scene_path,'s3lv_3D.txt');    
end

scene_raw_h = dlmread(fname_scene_h); xuse_h = scene_raw_h(:,1); yuse_h = scene_raw_h(:,2); zuse_h = scene_raw_h(:,3);
scene_raw_v = dlmread(fname_scene_v); xuse_v = scene_raw_v(:,1); yuse_v = scene_raw_v(:,2); zuse_v = scene_raw_v(:,3);

xuse = [xuse_h;xuse_v];
yuse = [yuse_h;yuse_v];
zuse = [zuse_h;zuse_v];

if CAM ~= 4
    ruse = sqrt(xuse.^2 + yuse.^2);
    puse = atan2(yuse,xuse);
    % [rs,phis,zs] = symmetrize_vec(ruse,puse,zuse,5,1);
    [rs,zs,phis] = symmetrize_and_move_to_hfp(ruse,zuse,puse,5,1,3,1);
    xuse = rs.*cos(phis);
    yuse = rs.*sin(phis);
    zuse = zs;
end
colors = lines;


for shot = myshots
    
    if CAM == 2
        file_endh = '_s2uh_1900_2000_ms.txt';
        file_endv = '_s2uv_1900_2000_ms.txt';
    elseif CAM == 3
        file_endh = '_s3lh_1900_2000_ms.txt';
        file_endv = '_s3lv_1900_2000_ms.txt';
    elseif CAM == 4
        file_endh = '_AEF40_s4lhdata.txt';
        file_endv = '_AEF40_s4lvdata.txt';
        % file_endh = '_s4lhdata_3000_3100.txt';
        % file_endv = '_s4lvdata_3000_3100.txt';
    end
    
    fnameir_h = fullfile(ir_path,strcat('20171129_0',num2str(shot),file_endh));
    fnameir_v = fullfile(ir_path,strcat('20171129_0',num2str(shot),file_endv));
    imageuse_h = dlmread(fnameir_h);
    imageuse_v = dlmread(fnameir_v);
    
    imageuse = [imageuse_h;imageuse_v];
    
    figure; hold on; box on;
    colorbar;
    title(['IR heat fluxshot ',num2str(shot)])
    
    ssize = ones(size(imageuse))*1;
    scatter3(xuse,yuse,zuse,ssize,imageuse,'.')
    iii = find(imageuse > 0.01);
    scatter3(xuse(iii),yuse(iii),zuse(iii),ssize(iii),imageuse(iii),'.')
    % Line checking
    %     [p1_blob,p2_blob,p1_edge,p2_edge] = define_W7X_divertor_1d_lines(0);
    %     [p1_blob,p2_blob,p1_edge,p2_edge,p1_blobd,p2_blobd,p1_blobd2,p2_blobd2,p1_blobd3,p2_blobd3]= define_W7X_divertor_1d_lines(0);
    [p1_edge,p2_edge,p1_edge2,p2_edge2,p1_edge3,p2_edge3,p1_blobd,p2_blobd,p1_blobd2,p2_blobd2,p1_blobd3,p2_blobd3] = define_W7X_divertor_1d_lines(0);
    blob_length = norm(p1_blobd-p2_blobd);
    edge_length = norm(p1_edge-p2_edge);
    icount_blob = 1;
    icount_edge = 1;
    dist_tol = 0.005;
    warning('should this be looped over xuse?')
    for i = 1:length(xuse_h)

        p0 = [xuse(i),yuse(i),zuse(i)];
        [dist,pu,CUTOFF,pu_unchecked,dist_unchecked,u] = distance_point_to_line_seg_3d(p1_blobd,p2_blobd,p0);
        if CUTOFF == 0 && dist < dist_tol
            psave_blob(icount_blob,:) = pu;
            dsave_blob(icount_blob) = imageuse(i);
            usave_blob(icount_blob) = u;
            icount_blob = icount_blob + 1;
        end
        [dist,pu,CUTOFF,pu_unchecked,dist_unchecked,u] = distance_point_to_line_seg_3d(p1_blobd2,p2_blobd2,p0);
        if CUTOFF == 0 && dist < dist_tol
            psave_blob(icount_blob,:) = pu;
            dsave_blob(icount_blob) = imageuse(i);
            usave_blob(icount_blob) = u;
            icount_blob = icount_blob + 1;
        end
        [dist,pu,CUTOFF,pu_unchecked,dist_unchecked,u] = distance_point_to_line_seg_3d(p1_blobd3,p2_blobd3,p0);
        if CUTOFF == 0 && dist < dist_tol
            psave_blob(icount_blob,:) = pu;
            dsave_blob(icount_blob) = imageuse(i);
            usave_blob(icount_blob) = u;
            icount_blob = icount_blob + 1;
        end
        [dist,pu,CUTOFF,pu_unchecked,dist_unchecked,u] = distance_point_to_line_seg_3d(p1_edge,p2_edge,p0);
        if CUTOFF == 0 && dist < dist_tol
            psave_edge(icount_edge,:) = pu;
            dsave_edge(icount_edge) = imageuse(i);
            usave_edge(icount_edge) = u;
            icount_edge = icount_edge + 1;
        end
        [dist,pu,CUTOFF,pu_unchecked,dist_unchecked,u] = distance_point_to_line_seg_3d(p1_edge2,p2_edge2,p0);
        if CUTOFF == 0 && dist < dist_tol
            psave_edge(icount_edge,:) = pu;
            dsave_edge(icount_edge) = imageuse(i);
            usave_edge(icount_edge) = u;
            icount_edge = icount_edge + 1;
        end
        [dist,pu,CUTOFF,pu_unchecked,dist_unchecked,u] = distance_point_to_line_seg_3d(p1_edge3,p2_edge3,p0);
        if CUTOFF == 0 && dist < dist_tol
            psave_edge(icount_edge,:) = pu;
            dsave_edge(icount_edge) = imageuse(i);
            usave_edge(icount_edge) = u;
            icount_edge = icount_edge + 1;
        end
        
    end
    
    
    
    [~,isort_edge] = sort(usave_edge);
    
    if shot == myshots(1)
        h1 = figure; hold on; box on;
        title(['Edge fluxshot ',num2str(shot)])
    else
        figure(h1);
    end
    if PLOT_UNAVERAGED
        plot(usave_edge(isort_edge)*edge_length,dsave_edge(isort_edge),'-','markersize',12)
    end
    
    NBIN = 100;
    x = usave_edge(isort_edge)*edge_length;
    y = dsave_edge(isort_edge);
    edges_want = linspace(0,max(x),NBIN);
    [nn,edges,bins] = histcounts(x,edges_want);
    clear YY
    for i = 1:length(nn)
        YY(i) = mean(y(bins == i));
    end
    plot((edges(1:end-1)+edges(2:end))/2,YY,'-','linewidth',3)
    
    % figure; hold on; box on;
    % plot(usave_edge(isort_edge),(dsave_edge(isort_edge)-dsave_edge(isort_edge(end)))./(max(dsave_edge)-max(dsave_edge(isort_edge(end)))),'kx-','markersize',12)
    % title('edge')
    
    [~,isort_blob] = sort(usave_blob);
    if shot == myshots(1)
        h2 = figure; hold on; box on;
        title(['blob fluxshot ',num2str(shot)])
    else
        figure(h2);
    end
    if PLOT_UNAVERAGED
        plot(usave_blob(isort_blob)*blob_length,dsave_blob(isort_blob),'-','markersize',12)
    end
    
    x = usave_blob(isort_blob)*blob_length;
    y = dsave_blob(isort_blob);
    edges_want = linspace(0,max(x),NBIN);
    [nn,edges,bins] = histcounts(x,edges_want);
    clear YY
    for i = 1:length(nn)
        YY(i) = mean(y(bins == i));
    end
    plot((edges(1:end-1)+edges(2:end))/2,YY,'-','linewidth',3)
    
    %
    % figure; hold on; box on;
    % plot(usave_blob(isort_blob),(dsave_blob(isort_blob)-dsave_blob(isort_blob(end-3)))./(max(dsave_blob)-max(dsave_blob(isort_blob(end-3)))),'kx-','markersize',12)
    % title('blob')
    
    
    
end
