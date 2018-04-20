clearvars;

% load('C:\Users\jjl\Dropbox (ORNL)\W7X\OP12a_Data\IR\shot27_irdata_cam10_t2_3_image1.mat')
% load('C:\Users\jjl\Dropbox (ORNL)\W7X\OP12a_Data\IR\shot27_irdata_cam11_t2_2.1_image1.mat')
% load('27_irdata_cam11_t2_2.1.mat')
ir_path = 'C:\Users\jjl\Dropbox (ORNL)\W7X\OP12a_Data\IR';
scene_path = 'C:\Users\jjl\Dropbox (ORNL)\W7X\scene_models';

PLOT_RAW = 0;

% myshots = [41,39,32,44,21];  % H 2MW density scan (0,11,22,32,43)
myshots = [10];
% myshots = [30,36];

% myshots = [40,37,27,43,10]; % He 2MW density scan (0,11,22,32,43)

% cam_strs = {'10','11','20','21','30','31','40','41','50','51'};
% cam_strs = {'11','20','40'};
cam_strs = {'40'}; Tmin = 350; Tmax = 650;
for shot = myshots
    for icam = 1:length(cam_strs)
        cam_str = cam_strs{icam};
        
        
        iframes = [1];
        %     cam_str = '10';  % 10,11
        time_str = 't3_3.1';
        fnameir = fullfile(ir_path,['shot',num2str(shot),'_irdata_cam',cam_str,'_',time_str,'.mat']);
        
        load(fnameir);
        if ir.status == 0
            disp(['Cam ',cam_str,' status 0 for shot ',num2str(shot)])
        else
            
            for iframe = iframes
                image = squeeze(ir.images(iframe,:,:)).';
                if PLOT_RAW
                    plot_w7x_IR_dynamic_clim(image.',Tmin,Tmax);
                    title(['Cam ',cam_str,' IR temperature shot ',num2str(shot)])
                end
                fname_scene = fullfile(scene_path,['AEF',cam_str,'\AEF',cam_str,'_scene_model.h5']);
%                 CAD = h5read(fname_scene,'/CAD');
                x = h5read(fname_scene,'/x'); y = h5read(fname_scene,'/y'); z = h5read(fname_scene,'/z');
                PFC = h5read(fname_scene,'/PFC');
                a = (PFC(:) == [1:15]) & image(:) > 350;
                iuse = any(a.');
                
                if 1
                    xuse = x(iuse);
                    yuse = y(iuse);
                    zuse = z(iuse);
                    imageuse = image(iuse);
                    figure; hold on; box on;
                    ssize = ones(size(imageuse))*0.5;
                    scatter3(xuse,yuse,zuse,ssize,imageuse,'.')
                    set(gca,'clim',[350,580])
                    colorbar;
                    title(['Cam ',cam_str,' IR temperature shot ',num2str(shot)])
                    
                    % Line checking
                    [p1_blob,p2_blob,p1_edge,p2_edge] = define_W7X_divertor_1d_lines();
                    icount_blob = 1;
                    icount_edge = 1;
                    dist_tol = 0.001;
                    for i = 1:length(xuse)
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
                end
%                 ruse = sqrt(xuse.^2 + yuse.^2);
%                 puse = atan2(yuse,xuse);
%                 nfp = 5;
%                 stellsym = 1;
%                 fp_want = 1;
%                 move_to_2nd_hfp = 0;
%                 [Rs1,Zs1,Ps1] = symmetrize_and_move_to_hfp(ruse,zuse,puse,nfp,stellsym,fp_want,move_to_2nd_hfp);
%                 Xs1 = Rs1.*cos(Ps1);
%                 Ys1 = Rs1.*sin(Ps1);
%                 
%                 figure; hold on; box on;
%                 scatter3(Xs1,Ys1,Zs1,[],iuse,'.')
%                 set(gca,'clim',[350,580])
%                 colorbar;
%                 
%                 afdadsaf
            end
        end
    end
    
    
end

[~,isort_edge] = sort(usave_edge);

figure; hold on; box on;
plot(usave_edge(isort_edge),dsave_edge(isort_edge),'kx-','markersize',12)
title('edge')

% figure; hold on; box on;
% plot(usave_edge(isort_edge),(dsave_edge(isort_edge)-dsave_edge(isort_edge(end)))./(max(dsave_edge)-max(dsave_edge(isort_edge(end)))),'kx-','markersize',12)
% title('edge')
% 
% [~,isort_blob] = sort(usave_blob);
% 
% figure; hold on; box on;
% plot(usave_blob(isort_blob),dsave_blob(isort_blob),'kx-','markersize',12)
% title('blob')
% 
% figure; hold on; box on;
% plot(usave_blob(isort_blob),(dsave_blob(isort_blob)-dsave_blob(isort_blob(end-3)))./(max(dsave_blob)-max(dsave_blob(isort_blob(end-3)))),'kx-','markersize',12)
% title('blob')

if 0
    CAD = h5read(fname_scene,'/CAD');
    x = h5read(fname_scene,'/x'); y = h5read(fname_scene,'/y'); z = h5read(fname_scene,'/z');
    PFC = h5read(fname_scene,'/PFC');
    % iuse = 1:numel(image);
    a = (PFC(:) == [1:15]) & image(:) > 350;
    iuse = any(a.');
    
    % image = fliplr(image);
    figure; hold on; box on;
    scatter3(x(iuse),y(iuse),z(iuse),[],image(iuse),'.')
    set(gca,'clim',[380,500])
    colorbar;
end