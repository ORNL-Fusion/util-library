clearvars;


run_path = 'C:\Work\DIII-D\164723\VMEC_XPAND\3059\';
fname = fullfile(run_path,'xpand_164723_3059.dat');

field = read_xpand_field_file(fname);

% R = 2.1;
% Z = -0.5;
% P = 10*pi/180;

% R = 2.1;
% P = 2*pi;
% Z = -1.5;

R = 2.1;
Z = 0.05;
P = 0.1;

tic;
interp3(field.z,field.r,field.phi,field.Br,Z,R,P)  % No idea why you have to swap these...

fprintf('This took %f\n',toc)
%%%%%%%%%%%%%%

iphi = floor(interp1(field.phi,1:field.nphi,P));
qq1 = interp2(field.r,field.z,field.Br(:,:,iphi).',R,Z);

phi_tmp = field.phi(field.nphi) - P;
if abs(phi_tmp) > eps
    qq2 = interp2(field.r,field.z,field.Br(:,:,iphi+1).',R,Z);
    interp1([field.phi(iphi),field.phi(iphi+1)],[qq1,qq2],P)
else
    qq1 
end


%%%%%%%%%%%%%%%%%
dr_grid = field.r(2) - field.r(1);
dz_grid = field.z(2) - field.z(1);
dphi_grid = field.phi(2) - field.phi(1);

ir = floor((R - field.r(1))/dr_grid) + 1;
iz = floor((Z - field.z(1))/dz_grid) + 1;

P = mod(P+2*pi,2*pi);
iphi = floor((P - field.phi(1))/dphi_grid) + 1;

% bi + lin
tic;
dr_grid = field.r(ir+1) - field.r(ir);
dz_grid = field.z(iz+1) - field.z(iz);
dphi_grid = field.phi(iphi+1) - field.phi(iphi);

Q1 = field.Br(ir:ir+1,iz:iz+1,iphi);
Q2 = field.Br(ir:ir+1,iz:iz+1,iphi+1);

drvec = [field.r(ir+1) - R,R - field.r(ir)];
dzvec = [field.z(iz+1) - Z;Z - field.z(iz)];
phi_fac = (P - field.phi(iphi))/dphi_grid;
(1-phi_fac)*drvec*Q1*dzvec/(dr_grid*dz_grid) + phi_fac*drvec*Q2*dzvec/(dr_grid*dz_grid)

fprintf('This took %f\n',toc)
