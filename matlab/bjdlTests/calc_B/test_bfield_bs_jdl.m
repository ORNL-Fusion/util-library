% Data
% clearvars;

% Shared variables
tolerance = 1e-9;
ntheta_per_wind = 100000;

% cur*sqrt(Brg^2 + Bzg^2);
% Brg = cur*Brg;
%% Test singular circular coil (Analytic)
r1 = 0.1;
r2 = 0.11;
z1 = 0.5;
dz = 0.1;
nturns = 1;
nlayers = 1;
cur = 1e3;



[coil,current] = build_circular_coil(r1,r2,z1,dz,nturns,nlayers,cur,ntheta_per_wind);

P_x = 0.09;
P_y = 0.01;
P_z = z1+dz;
P_r = sqrt(P_x^2 + P_y^2);

[Bx,By,Bz,Btot]=bfield_bs_jdl(P_x,P_y,P_z,coil,current);
Br = sqrt(Bx^2+By^2);


fw = dz/nturns; % Width (axial length) of each turn
fh = (r2-r1)/nlayers; % Thickness (radial length) of each layer

% r,z positions of the center of each filament
zwind0 = repmat(linspace(z1+fw/2,z1+dz-fw/2,nturns),1,nlayers);
rwind0 = repmat(linspace(r1+fh/2,r2-fh/2,nlayers),1,nturns);
[Brg,Bzg]=bfield_circular_coil_analytic(rwind0,zwind0,P_r,P_z);

errorTest = abs(cur*Brg - Br);
assert(errorTest < tolerance,'Did not match Analytic circular coil result. Err = %.2e',errorTest)


%% Test stacked circular coil (Analytic)
 
r1 = 0.1;
r2 = 0.11;
z1 = 0.5;
dz = 0.1;
nturns = 3;
nlayers = 4;
cur = 1e3;

[coil,current] = build_circular_coil(r1,r2,z1,dz,nturns,nlayers,cur,ntheta_per_wind);

P_x = 0.09;
P_y = 0.01;
P_z = z1+dz;
P_r = sqrt(P_x^2 + P_y^2);

fw = dz/nturns; % Width (axial length) of each turn
fh = (r2-r1)/nlayers; % Thickness (radial length) of each layer

% r,z positions of the center of each filament
zwind0 = repmat(linspace(z1+fw/2,z1+dz-fw/2,nturns),1,nlayers);
rwind0 = repmat(linspace(r1+fh/2,r2-fh/2,nlayers),1,nturns);

[Bx,By,Bz,Btot]=bfield_bs_jdl(P_x,P_y,P_z,coil,current);
Br = sqrt(Bx^2+By^2);


for i = 1:numel(rwind0)
    [Brg(i),Bzg(i)]=bfield_circular_coil_analytic(rwind0(i),zwind0(i),P_r,P_z);
end

errorTest = abs(cur*sum(Brg) - Br);
assert(errorTest < tolerance,'Did not match Analytic circular coil result. Err = %.2e',errorTest)

