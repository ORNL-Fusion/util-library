clearvars;

debug_plots = 1;

mypath = fileparts(mfilename('fullpath'));  
fname_save = fullfile(mypath,'fs_data.mat');

ves = load_W7X_vessel(0,0);
figure; hold on; box on;
for i = 1:5
    vx = ves.r.*cos(ves.phi+(i-1)*2*pi/5);
    vy = ves.r.*sin(ves.phi+(i-1)*2*pi/5);
    vz = ves.z;
    plot3(vx,vy,vz,'b')
    plot3(vx.',vy.',vz.','b')
end
fs = define_fs_raw(0,0);

names = fieldnames(fs);
nfs = length(names);
L_extend = 4;
for i = 1:nfs
    ifs = getfield(fs,names{i});
    v1 = ifs.p2 - ifs.p1;
    L1 = sqrt(sum(v1.^2));
    n1 = v1./L1;
    p22 = ifs.p1 + L_extend*n1;
    plot3(ifs.p1(1),ifs.p1(2),ifs.p1(3),'ko')
    plot3(ifs.p2(1),ifs.p2(2),ifs.p2(3),'ro')
    plot3([ifs.p1(1),p22(1)],[ifs.p1(2),p22(2)],[ifs.p1(3),p22(3)],'m')
    text(ifs.p1(1),ifs.p1(2),ifs.p1(3),names{i})
end

for i = 1:nfs
    ifs = getfield(fs,names{i});
    v1 = ifs.p2 - ifs.p1;
    L1 = sqrt(sum(v1.^2));
    n1 = v1./L1;
    p22 = ifs.p1 + L_extend*n1;

    dL_fine = 0.001;    
    Lmax = 4;
%     nL = round(Lmax/dL);
%     Ltest = linspace(0,Lmax,Lmax/dL);
    
    % now find vessel intersections
    if abs(atan2(ifs.p1(2),ifs.p1(1)) - atan2(ifs.p2(2),ifs.p2(1))) < 1e-6
        const_phi = 1;
    else
        const_phi = 0;
    end
    fprintf(strcat(names{i},', const_phi = ',num2str(const_phi)))
    fprintf('.  phi = [%f, %f] [deg]\n',atan2(ifs.p1(2),ifs.p1(1))*180/pi,atan2(ifs.p2(2),ifs.p2(1))*180/pi)
    fs = setfield(fs,names{i},'const_phi',const_phi);
    if debug_plots
        figure; hold on; box on;
        title(strcat(names{i},', const\_phi = ',num2str(const_phi)))
    end
    change_count = 0;
    isin_last = 0;
    Ltest = 0;
    while Ltest <= Lmax
        p22 = ifs.p1 + Ltest*n1;
        xtest = p22(1);
        ytest = p22(2);
        ztest = p22(3);
        rtest = sqrt(xtest^2 + ytest^2);
        ptest = atan2(ytest,xtest);
        if ~const_phi || Ltest == 0
            ves_cut = cut_W7X_vessel(ves,ptest);
        end
        isin = inpolygon(rtest,ztest,ves_cut.r,ves_cut.z);
        if debug_plots            
            if isin == 1
                plot(rtest,ztest,'ro')
                plot(ves_cut.r,ves_cut.z,'b')
            else
                plot(rtest,ztest,'kx')
                plot(ves_cut.r,ves_cut.z,'k')
            end
        end
                
        if isin ~= isin_last
            change_count = change_count + 1;
            if change_count == 1                
                fs = setfield(fs,names{i},'p1_ves',[xtest,ytest,ztest]);
                fs = setfield(fs,names{i},'p1cyl_ves',[rtest,ptest,ztest]);
            end
            if change_count > 1
                fs = setfield(fs,names{i},'p2_ves',[xtest,ytest,ztest]);
                fs = setfield(fs,names{i},'p2cyl_ves',[rtest,ptest,ztest]);
                break;
            end
            isin_last = isin;
        end
        min_dL = min(sqrt((ves_cut.r-rtest).^2 + (ves_cut.z-ztest).^2));
        if min_dL > dL_fine
            dL = min_dL/2;
        else
            dL=dL_fine;
        end
        Ltest = Ltest + dL;
    end
end
    
save(fname_save,'fs')


plotit = 1;
if plotit
    ves = load_W7X_vessel;
    names = fieldnames(fs);
    nfs = length(names);    
    for i = 1:nfs
        ifs = getfield(fs,names{i});
        figure; hold on; box on;
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
    end
end