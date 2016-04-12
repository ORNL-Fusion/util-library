function fs = get_filterscope_geo

mypath = fileparts(mfilename('fullpath'));  
fname = fullfile(mypath,'fs_data.mat');
load(fname);

plotit = 0;
subit = 0;
if plotit
    ves = load_W7X_vessel;
    names = fieldnames(fs);
    nfs = length(names);        
    for i = 1:nfs
        ifs = getfield(fs,names{i});        
        if i == 1 || ~subit
            figure; hold on; box on;
        end
        if subit
            subplot(3,2,i); hold on; box on;
        end
        phi1 = atan2(ifs.p1_ves(2),ifs.p1_ves(1));
        r1 = sqrt(ifs.p1_ves(1)^2 + ifs.p1_ves(2)^2);
        z1 = ifs.p1_ves(3);
        phi2 = atan2(ifs.p2_ves(2),ifs.p2_ves(1));
        r2 = sqrt(ifs.p2_ves(1)^2 + ifs.p2_ves(2)^2);
        z2 = ifs.p2_ves(3);
        
        ves_cut = cut_W7X_vessel(ves,(phi1+phi2)/2);
        
        plot(ves_cut.r,ves_cut.z,'b')
        plot([r1,r2],[z1,z2],'r')
%         plot([ifs.p1cyl_ves(1),ifs.p2cyl_ves(1)],[ifs.p1cyl_ves(3),ifs.p2cyl_ves(3)],'r')
        title(sprintf('%s\n \\phi = [%4.1f,%4.1f]',names{i},phi1*180/pi,phi2*180/pi))        
        axis equal;
        xlabel('R [m]')
        ylabel('Z [m]')
    end
end