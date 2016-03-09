clearvars;
[coil,COIL] = make_larry_coils(100);

r_start = -0.07;
r_end = 0.07;
z_start = -0;
z_end = 5;

nr = 100;
nz = 1000;

Rstart = 0.05;
Zstart = 0;
phistart = 0;
dl = 0.005;
nsteps = 1000;
nowarn = 0;
type = 'just_coils';
rmp = COIL;
g = [];
[s,ierr,i_last_good]=follow_fieldlines_rzphi_dl(g,rmp,Rstart,Zstart,phistart,dl,nsteps,nowarn,'just_coils');

figure; hold on; box on;
plot(s.r,s.z)

adsfdsaf

r_array = linspace(r_start,r_end,nr);
z_array = linspace(z_start,z_end,nz);

[rg,zg] = meshgrid(r_array,z_array);
for ir = 1:nr
    fprintf('Working on ir %i of %i\n',[ir,nr]);
    for iz = 1:nz
        x = r_array(ir);
        z = z_array(iz);
        y = 0;
        
        [Bx(ir,iz),By(ir,iz),Bz(ir,iz)]=bfield_bs_jdl(x,y,z,COIL.coil,COIL.current);        
    end
end

modB = sqrt(Bx.^2 + By.^2 + Bz.^2);


if 0
    figure; hold on; box on;
    surf(rg.',zg.',modB,'edgecolor','none')
    
    
    figure; hold on; box on;
    plot(z_array,modB(end/2,:))
end

% r_array = linspace(r_start,r_end,nr);
% z_array = linspace(z_start,z_end,nz);
% 
% [rg,zg] = meshgrid(r_array,z_array);
% for ir = 1:nr
%     fprintf('Working on ir %i of %i\n',[ir,nr]);
%     for iz = 1:nz
%         r = r_array(ir);
%         z = z_array(iz);
%         
%         [Br(ir,iz),Bphi(ir,iz),Bz(ir,iz)]=bfield_bs_cyl(r,0,z,coil.coil,coil.current);        
%     end
% end
% 
% modB = sqrt(Br.^2 + Bphi.^2 + Bz.^2);
% 
% figure; hold on; box on;
% contour(rg.',zg.',modB)
% 
% % Rstart = 0.05;
% % Zstart = 0;
% % phistart = 0;
% % dphi = 0.1*pi/180;
% % nsteps = 100;
% % nowarn = 0;
% % type = 'just_coils';
% % g = [];
% % [s,ierr,i_last_good]=follow_fieldlines_rzphi(g,coil,Rstart,Zstart,phistart,dphi,nsteps,nowarn,type);
% % 
% % x = s.r.*cos(s.phi);
% % y = s.r.*sin(s.phi);
% % z = s.z;
% % 
% % plot3(x,y,z)