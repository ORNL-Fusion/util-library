clearvars;

ir_path = 'C:\Users\jjl\Dropbox (ORNL)\W7X\OP12a_Data\Jeremy_H';
scene_path = ir_path;

% myshots = [41,39,32,44,21];  % H 2MW density scan (0,11,22,32,43)
myshots = [32];

fname_scene_h = fullfile(scene_path,'s4lh_3Dlines_coordinate.txt');
scene_raw_h = dlmread(fname_scene_h); xuse_h = scene_raw_h(:,1); yuse_h = scene_raw_h(:,2); zuse_h = scene_raw_h(:,3);
fname_scene_v = fullfile(scene_path,'s4lv_3Dlines_coordinate.txt');
scene_raw_v = dlmread(fname_scene_v); xuse_v = scene_raw_v(:,1); yuse_v = scene_raw_v(:,2); zuse_v = scene_raw_v(:,3);

xuse = [xuse_h;xuse_v];
yuse = [yuse_h;yuse_v];
zuse = [zuse_h;zuse_v];

colors = lines;


for shot = myshots
    
    fnameir_h = fullfile(ir_path,strcat('20171129_0',num2str(shot),'_AEF40_s4lhdata.txt'));
    imageuse_h = dlmread(fnameir_h);
    fnameir_v = fullfile(ir_path,strcat('20171129_0',num2str(shot),'_AEF40_s4lvdata.txt'));
    imageuse_v = dlmread(fnameir_v);
    
    imageuse = [imageuse_h;imageuse_v];
    
    figure; hold on; box on;
    colorbar;
    title(['IR heat fluxshot ',num2str(shot)])
    
    ssize = ones(size(imageuse))*1;
    scatter3(xuse,yuse,zuse,ssize,imageuse,'.')       
%     iii = find(imageuse > 0.01);
%     scatter3(xuse(iii),yuse(iii),zuse(iii),ssize(iii),imageuse(iii),'.')
    % Line checking
    [p1_blob,p2_blob,p1_edge,p2_edge] = define_W7X_divertor_1d_lines(0);
    icount_blob = 1;
    icount_edge = 1;
    dist_tol = 0.005;
    for i = 1:length(xuse_h)
        p0 = [xuse(i),yuse(i),zuse(i)];
        [dist,pu,CUTOFF,pu_unchecked,dist_unchecked,u] = distance_point_to_line_seg_3d(p1_blob,p2_blob,p0);
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
        
    end
    
    
    
    [~,isort_edge] = sort(usave_edge);
    
    if shot == myshots(1)
        h1 = figure; hold on; box on;
        title(['Edge fluxshot ',num2str(shot)])
    else
        figure(h1);
    end
    
    plot(usave_edge(isort_edge),dsave_edge(isort_edge),'-','markersize',12)
    
    %
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
    plot(usave_blob(isort_blob),dsave_blob(isort_blob),'-','markersize',12)
    
    %
    % figure; hold on; box on;
    % plot(usave_blob(isort_blob),(dsave_blob(isort_blob)-dsave_blob(isort_blob(end-3)))./(max(dsave_blob)-max(dsave_blob(isort_blob(end-3)))),'kx-','markersize',12)
    % title('blob')
    
    
    
end
