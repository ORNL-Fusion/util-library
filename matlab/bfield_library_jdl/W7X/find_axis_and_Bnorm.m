clearvars;


coils_file = 'C:\Work\Stellarator\W7X EMC3 modeling\Mark coils and input\coils.w7x';
coil = load_vmec_coils_file(coils_file);
winding_array = [108,108,108,108,108,36,36,8,8];

% taper_norm = [12130, 12000, 12255, 13541, 13635, 9000, -2900, 0, 0].*winding_array; Inorm = 1;  % 22kA OP2
taper_norm = [1, 1.02, 1.08, 0.97, 0.88, 0.15, -0.15, 0, 0].*winding_array; Inorm = 12556;  % Narrow mirror
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.2200, -0.0800, 0.0150, -0.0150]; Inorm = 1340931;%  0kA mimic
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1950, -0.1050, 0.0150, -0.0150]; Inorm = 1353630; % 11kA mimic
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1700, -0.1300, 0.0150, -0.0150]; Inorm = 1366546; % 22kA mimic  -- corrected IS1
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1450, -0.1550, 0.0150, -0.0150]; Inorm = 1379685; % 32kA mimic
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1200, -0.1800, 0.0150, -0.0150]; Inorm = 1393054; % 43kA mimic
taper = Inorm*taper_norm./winding_array;

coil = set_w7x_current(coil,taper); % taper = [I1,I2,I3,I4,I5,IA,IB,IS1,IS2];
bfield.type = 'just_coils';
bfield.coil = coil.coil;
bfield.current = coil.current;
bfield.nsym = coil.num_periods;

phistart = 0*pi/180; Rax_guess = 5.92; Zax_guess = 0;
% phistart = 36*pi/180; Rax_guess = 5.16; Zax_guess = 0;


[Rax,Zax] = find_axis(bfield,phistart,Rax_guess,Zax_guess);
[Bax] = bfield_general_rzphi(Rax,Zax,phistart,bfield);
B0 = sqrt(Bax.br.^2 + Bax.bz.^2 + Bax.bphi.^2);
fprintf('B0 = %f\n',B0)



