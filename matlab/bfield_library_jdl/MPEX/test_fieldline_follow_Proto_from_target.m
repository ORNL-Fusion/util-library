% function test_fieldline_follow_Proto
clearvars;
x0_guess = []; y0_guess = []; force_guess = 0;
verbose = 0;
if 0
    
    helicon_current = 500;
    current_A = 6400;
    current_B = 0;
    current_C = 1;
    config = 'flat';
    skimmer = 1;
    shot = 0;
    plasma_radius_cm = 1;
    target_position = 2;
    sleeve = 1;
    
else
    
%     shot = 7418;  mytitle = 'I_A = 6368 A, no skimmer';  x0_guess = -.5512;y0_guess = -2.533; force_guess = 1; %shots = 7400 + [0,3:6,8,10,12:13,16,17,18];
%      shot = 7477; mytitle = 'I_A = 3300 A, with skimmer';x0_guess = -.5877;y0_guess = -2.8914; force_guess = 1; % shots = 7400 + [77,87,88,92:98,100,101,103];
    shot = 7277;  mytitle = 'I_A = 3300 A, no skimmer';  x0_guess = -0.604712; y0_guess = -3.026471; force_guess = 0; %shots = 7400 + [0,3:6,8,10,12:13,16,17,18];
%     shot = 7278;  mytitle = 'I_A = 3300 A, no skimmer';  x0_guess = -0.604712; y0_guess = -3.026471; force_guess = 0; %shots = 7400 + [0,3:6,8,10,12:13,16,17,18];
%     shot = 7674;  mytitle = 'I_A = 3300 A, no skimmer';  x0_guess = -0.591378; y0_guess = -2.477497; force_guess = 0; %shots = 7400 + [0,3:6,8,10,12:13,16,17,18];
    
%     shot = 6547; x0_guess = 0.05; y0_guess = -1.7; force_guess = 0;
    if isempty(x0_guess)
        [rr_cm_IR,dd_cm_IR,plasma_radius_cm] = fit_IR_data(shot,1,0,0);
    else
        [rr_cm_IR,dd_cm_IR,plasma_radius_cm] = fit_IR_data(shot,1,x0_guess,y0_guess,force_guess);
    end
    [helicon_current,current_A,current_B,config,skimmer,current_C] = get_Proto_current(shot);
    target_position = 1;
    sleeve = 0;
    
end
drawnow;

[coil,current] = build_Proto_coils(helicon_current,current_A,current_B,config,verbose,current_C);
geo = get_Proto_geometry(0,0,skimmer,target_position,sleeve);

bfield.coil = coil;
bfield.current = current;
bfield.type = 'just_coils';
bfield.vessel_clip_r   = geo.vessel_clip_r;
bfield.vessel_clip_z   = geo.vessel_clip_z;
bfield.vessel_clip_phi = 0;
bfield.stop_at_vessel = 1;

num_lines = 25; 
% rr = linspace(1e-3,0.04,num_lines);
rr = linspace(1e-3,2.0*plasma_radius_cm/100,num_lines);  % Initial R positions
zz = geo.target.z*ones(size(rr));                        % Initial Z positions
L = 2.5;                                                 % Length to follow field lines
dl = -0.01;                                              % Step size (m)
nsteps = abs(L/dl);
phistart = zeros(size(rr));
tic;
f2 = follow_fieldlines_rzphi_dz(bfield,rr,zz(1),phistart,dl,nsteps);
for i =1:length(zz)
    f.r = f2.r(:,i); f.z = f2.z(:,i); f.phi = f2.phi(:,i);
    fsave{i} = f;
end
toc
% tic;
% for i = 1:length(zz)
%     fprintf('Line %d of %d\n',i,num_lines)
%     f = follow_fieldlines_rzphi_dz(bfield,rr(i),zz(i),phistart(i),dl,nsteps);
%     plot(f.z,f.r,'b','linewidth',2)
%     fsave{i} = f;
% end
% toc

figure; hold on; box on;
for i =1:length(zz)
    plot(fsave{i}.z,fsave{i}.r,'b','linewidth',2)
end
set(gca,'fontsize',14)
xlabel('Z [m]','fontsize',14)
ylabel('R [m]','fontsize',14)
title(['Shot ',num2str(shot)])
geo = get_Proto_geometry(1,0,skimmer,target_position,sleeve);
axis([0.5,3.5,0,0.15])

figure; hold on; box on;
xlabel('Z [m]','fontsize',14)
ylabel('R [m]','fontsize',14)
set(gca,'fontsize',14)
title(['Shot ',num2str(shot)])
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

    IR_interp = interp1(rr_cm_IR/100,dd_cm_IR,reval);
    patch(z,r,IR_interp,'edgecolor','none')
end
colorbar;
[rcoil,zcoil] = get_coil_cross_sections;
for i = 1:size(rcoil,1)
    plot(zcoil(i,:),rcoil(i,:),'r')
end
% geo = get_Proto_geometry(1,0,skimmer);
geo = get_Proto_geometry(1,0,skimmer,target_position,sleeve);
plot(geo.target.z*[1,1],geo.target.r*[0,1],'k','linewidth',3)
plot([geo.helicon.z1,geo.helicon.z2],geo.helicon.r*[1,1],'k','linewidth',3)
axis([0.5,3.5,0,0.15])
title(sprintf('Shot %d, I_H^*=%3.0f A',shot,helicon_current*(3300/current_A)))


% CLIP AT VESSEL
figure; hold on; box on;
xlabel('Z [m]','fontsize',14)
ylabel('R [m]','fontsize',14)
set(gca,'fontsize',14)
geo = get_Proto_geometry(1,0,skimmer,target_position,sleeve);
colorbar;
title(['Shot ',num2str(shot)])


for i = 1:num_lines    
    if i == 1        
        z = [fsave{i}.z;flipud(fsave{i}.z)];
        r = [fsave{i}.r;0*flipud(fsave{i}.r)]; 
        reval = (r(1)+r(end))/2;
    else
        
        rout = fsave{i}.r;
        zout = fsave{i}.z;        
        isin2 = inpolygon(rout,zout,geo.vessel_clip_r,geo.vessel_clip_z);
        is2 = find(isin2 == 0,1,'first')-1;
        if isempty(is2)
            is2 = length(isin2);
        end                  
        rout = rout(1:is2);
        zout = zout(1:is2);
        Lout=sum(sqrt(diff(rout).^2+diff(zout).^2));    
        
        
        rin = fsave{i-1}.r;
        zin = fsave{i-1}.z;        
        isin1 = inpolygon(rin,zin,geo.vessel_clip_r,geo.vessel_clip_z);
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
     
    IR_interp = interp1(rr_cm_IR/100,dd_cm_IR,reval);
    patch(z,r,IR_interp,'edgecolor','none')
end
axis([0.5,3.5,0,0.15])
title(sprintf('Shot %d, I_H^*=%3.0f A',shot,helicon_current*(3300/current_A)))
