clearvars;
run_path = 'C:\Work\FLARE\2dts_scan_as_4000ms\';

nz = 3;

for iz = 1:nz
    data = dlmread([run_path,'base_grid_',num2str(iz-1),'.dat'],'',4,0);
    rgrid{iz}(:) = data(:,1);
    zgrid{iz}(:) = data(:,2);
end



figure; hold on; box on;
colormap(colorflipper);
cf = colorflipper(4);
for iz = 1:nz
    plot(rgrid{iz},zgrid{iz},'.','color',cf(iz,:))
end

