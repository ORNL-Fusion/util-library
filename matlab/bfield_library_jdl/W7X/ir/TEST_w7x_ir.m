clearvars;

% load('C:\Users\jjl\Dropbox (ORNL)\W7X\OP12a_Data\IR\shot27_irdata_cam10_t2_3_image1.mat')
load('C:\Users\jjl\Dropbox (ORNL)\W7X\OP12a_Data\IR\shot27_irdata_cam11_t2_2.1_image1.mat')

image = ir.images.';
plot_w7x_IR_dynamic_clim(image,350,500);
title('IR temperature')

f = 'C:\Users\jjl\Dropbox (ORNL)\W7X\scene_models\AEF11\AEF11_scene_model.h5';  
% f = 'C:\Users\jjl\Dropbox (ORNL)\W7X\scene_models\AEF10\AEF10_scene_model.h5';MAKE SURE IT MATCHES IMAGE!!!!
% f2 = 'C:\Users\jjl\Dropbox (ORNL)\W7X\scene_models\AEF10\testimage.h5';

% image = squeeze(h5read(f2,'/images'));

% h5disp(f)
CAD = h5read(f,'/CAD');
x = h5read(f,'/x'); y = h5read(f,'/y'); z = h5read(f,'/z'); 
PFC = h5read(f,'/PFC');
% 
% plot_w7x_IR_dynamic_clim(CAD);
% title('CAD')

% plot_w7x_IR_dynamic_clim(PFC);

% plot3(x,y,z,'.')

% iuse = 1:numel(image);
a = (PFC(:) == [1:15]) & image(:) > 400;
iuse = any(a.');

% image = fliplr(image);
figure; hold on; box on;
scatter3(x(iuse),y(iuse),z(iuse),[],image(iuse),'.')
set(gca,'clim',[350,500])