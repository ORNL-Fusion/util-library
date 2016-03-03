% function test_fieldline_follow_Proto
clearvars;
% helicon_current = -80;
% current_A = 6607;
% current_B = 6607;
% config = 'focus';

helicon_current = -70;
current_A = 3300;
current_B = 0;
config = 'flat';

[coil,current] = build_Proto_coils(helicon_current,current_A,current_B,config);

bfield.coil = coil;
bfield.current = current;
bfield.type = 'just_coils';

Rstart = 0.01;
Zstart = 0;
phistart = 0;
dl = -0.01;
L = 5;
nsteps = abs(5/dl);

% f = follow_fieldlines_rzphi_dl(bfield,Rstart,Zstart,phistart,dl,nsteps);

figure; hold on; box on;
% plot(f.z,f.r,'b','linewidth',2)
xlabel('Z [m]','fontsize',14)
ylabel('R [m]','fontsize',14)
set(gca,'fontsize',14)
axis([0,5,0,0.2])

[rcoil,zcoil] = get_coil_cross_sections;
for i = 1:size(rcoil,1)
    plot(zcoil(i,:),rcoil(i,:),'r')
end


geo = get_Proto_geometry;
plot(geo.target.z*[1,1],geo.target.r*[0,1],'k','linewidth',3)
plot([geo.helicon.z1,geo.helicon.z2],geo.helicon.r*[1,1],'k','linewidth',3)

num_lines = 10;
% rr = linspace(1e-3,0.04,num_lines);
rr = linspace(1e-3,0.01,num_lines);
zz = geo.target.z*ones(size(rr));
L = 3;
nsteps = abs(L/dl);
for i = 1:length(zz)
    fprintf('Line %d of %d\n',i,num_lines)
    f = follow_fieldlines_rzphi_dz(bfield,rr(i),zz(i),phistart,dl,nsteps);
    plot(f.z,f.r,'b','linewidth',2)
    fsave{i} = f;
end



[rr_cm,dd] = plot_IR_data_raw(6547);

for i = 1:num_lines    
    if i == 1
        z = [fsave{i}.z;flipud(fsave{i}.z)];
        r = [fsave{i}.r;0*flipud(fsave{i}.r)]; 
        reval = (r(1)+r(end))/2;
    else
        z = [fsave{i-1}.z;flipud(fsave{i}.z)];
        r = [fsave{i-1}.r;flipud(fsave{i}.r)];
        reval = (r(1)+r(end))/2;
    end

    IR_interp = interp1(rr_cm/100,dd,reval);
    patch(z,r,IR_interp,'edgecolor','none')
end
geo = get_Proto_geometry(1,0);
colorbar;




% CLIP AT VESSEL
vessel_clip_z = [0,geo.vessel.z,5];
vessel_clip_r = [0,geo.vessel.r,0];

figure; hold on; box on;
xlabel('Z [m]','fontsize',14)
ylabel('R [m]','fontsize',14)
set(gca,'fontsize',14)
axis([0,5,0,0.2])
geo = get_Proto_geometry(1,0);
colorbar;


for i = 1:num_lines    
    if i == 1        
        z = [fsave{i}.z;flipud(fsave{i}.z)];
        r = [fsave{i}.r;0*flipud(fsave{i}.r)]; 
        reval = (r(1)+r(end))/2;
    else
        
        rout = fsave{i}.r;
        zout = fsave{i}.z;        
        isin2 = inpolygon(rout,zout,vessel_clip_r,vessel_clip_z);
        is2 = find(isin2 == 0,1,'first')-1;
        if isempty(is2)
            is2 = length(isin2);
        end                  
        rout = rout(1:is2);
        zout = zout(1:is2);
        Lout=sum(sqrt(diff(rout).^2+diff(zout).^2));    
        
        
        rin = fsave{i-1}.r;
        zin = fsave{i-1}.z;        
        isin1 = inpolygon(rin,zin,vessel_clip_r,vessel_clip_z);
        is1 = find(isin1 == 0,1,'first')-1;
        if isempty(is1)
            is1 = length(isin1);
        end        
        rin = rin(1:is1);
        zin = zin(1:is1);

        Lin=sum(sqrt(diff(rin).^2+diff(zin).^2));
        if Lin > Lout
            [icurve_near_L,err_near_L,R_L,Z_L] = move_L_on_C(Lout,rin,zin);
            rin = rin(1:icurve_near_L);
            zin = zin(1:icurve_near_L);
        end                
        z = [zin;flipud(zout)];
        r = [rin;flipud(rout)];
        reval = (r(1)+r(end))/2;
        
        
    end
     
    IR_interp = interp1(rr_cm/100,dd,reval);
    patch(z,r,IR_interp,'edgecolor','none')
end

