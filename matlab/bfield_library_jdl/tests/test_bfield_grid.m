clearvars;



temp = load('C:\Work\Stellarator\W7X EMC3 modeling\Mark coils and input\Bgrid_w7x_0kA_mimic_180x164x132.mat');
Bgrid = temp.Bgrid;
nowarn = 0;

coils_file = 'C:\Work\Stellarator\W7X EMC3 modeling\Mark coils and input\coils.w7x';
coil = load_vmec_coils_file(coils_file);
winding_array = [108,108,108,108,108,36,36,8,8];
taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.2200, -0.0800, -0.0150, 0.0150]; Inorm = 1.341e6;  %  0kA mimic
taper = Inorm*taper_norm./winding_array;
coil = set_w7x_current(coil,taper); % taper = [I1,I2,I3,I4,I5,IA,IB,IS1,IS2];
bfield.type = 'just_coils';
bfield.coil = coil.coil;
bfield.current = coil.current;
bfield.nsym = coil.num_periods;


Rmin = 4.25;
Rmax = 6.45;
Zmin = -1.2;
Zmax = 1.2;


% num_test = 20000;
num_test = 1;

t_grid = 0;
t_BS = 0;
for i = 1:num_test
    r1 = rand([3,1]);
    R = (Rmax-Rmin)*r1(1)+Rmin;
    Z = (Zmax-Zmin)*r1(2)+Zmin;
    P_rad = 2*pi*r1(3);
%    R = 5.186540; 
%    Z = -0.605598;
%    P_rad = 298.546024*pi/180;
    tic;
    [Br,Bz,Bphi]=bfield_grid(R,Z,P_rad,Bgrid,nowarn);
    t_grid = t_grid + toc;
    
    % [Bout,ierr] = bfield_general_rzphi(R,Z,P_rad,bfield,nowarn);
    
    tic;
    [Br_BS,Bphi_BS,Bz_BS]=bfield_bs_cyl(R,P_rad,Z,bfield.coil,bfield.current,nowarn);
    t_BS = t_BS + toc;
%     
    if num_test == 1
        fprintf('B from grid, B-S\n')
        fprintf('Evaluated at [R,Z,phi (deg)] = [%f, %f, %f]\n',R,Z,P_rad*180/pi)
        fprintf('Br   = %f, %f \n',Br,Br_BS)
        fprintf('Bz   = %f, %f \n',Bz,Bz_BS)
        fprintf('Bphi = %f, %f \n',Bphi,Bphi_BS)
        fprintf('Time for grid calc: %f\n',t_grid)
        fprintf('Time for B-S calc : %f\n',t_BS)
    end
end
fprintf('Time for grid calc: %f\n',t_grid)
fprintf('Time for B-S calc : %f\n',t_BS)
fprintf('Speedup? : %f\n',t_BS/t_grid)
% Now do splines?
