clearvars;

%% Definitions
% Throughput Q 
% (Torr-Liters/sec) = P (Torr) * S (L/s) 
% (Pa m^3/s)
% also called "standard condition volumetric flow"
%
% 

% 0C = 273.15 K

pc = phys_const;

% Ideal gas law
% p*V = (m/M)*R*T = (N/NA)*R*T = N*k*T
% p = n*k*T
% m = mass
% M = molar mass (kg/mol)
% N = number of particles
% n = density
% NA: 1/mol
% R = NA*k: J/K/mol General gas constant
% kB: J/K

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Temperature

    % Define in deg C
%     T_C = 20;
%     T_K = T_C + 273.15;

%     % Define in K
    T_K = 273.15;  % 0C
    T_C = T_K - pc.convert.C_to_K_Add;
    T_eV = T_K*pc.convert.K_to_eV;
    
%     % Define in eV
%     T_eV = 0.025;
%     T_K = T_eV*pc.e/pc.kB


%% Mass

    % Define as molar mass
%     M_gPermol = 2.014;
    M_gPermol = 4.002602;
    M_kgPermol = M_gPermol/1e3;
    m = M_kgPermol/pc.NAvogadro;

%% Pressure

    % Define in Pa
    p_Pa = 1e3*100; 
    p_mTorr = p_Pa*pc.convert.Pa_to_mTorr;

% Molar masses (g/mol)
% D2 = 2.014
% H2 = 2.016
% He = 4.002602
% Water vapor = 18.02
% Neon = 20.179
% N2 = 28.0134
% Air = 28.966
% Ar = 39.948
% CO2 = 44.01

% Radii in A
% H0 - 0.529
% H2 - 1.37   (I've also seen 1.2)
% N2 - 1.75
% He - 1.1


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cbar  = sqrt(8*pc.kB*T_K/pi/m);
cprob = sqrt(2*pc.kB*T_K/m);

n = p_Pa*pc.NAvogadro/(pc.Rgas*T_K); % should be the same as p/(kT)
rho = m*n; 

fprintf('For gas with: \n')
fprintf('             Molar mass = %.3f (g/mol), m = %.3e (kg)\n',M_gPermol,m)
fprintf('             T = %.3f (K) = %.3f (C), = %.3e (eV)\n',T_K,T_C,T_eV)
fprintf('             p = %.3f (Pa) = %.3f (mTorr)\n',p_Pa,p_mTorr)
fprintf('  cbar = %.3e (m/s)\n',cbar)
fprintf('  n = %.3e (m^-3), rho = %.3e (kg/m^3)\n',n,rho)


% Given a flow rate of molecular hydrogen, assume the pressure
% is for molecular particles at 273.15 K. Then divide by kT and appropriate 
% pV unit conversion to get particles/s, where particles = molecules. Then if you
% want atoms/s multiply this by 2. So convert Torr-L/s to Pa*m^3/s by dividing by
% 7.5.  Then divide by kT, this gives 1 Torr-L/s = 3.5355e19 molecules per second 
% and 7.07e19 atoms/s, and 11.33 A. 


% r_m = 1.37e-10;  % radius 

%
% d_m = 2*r_m;
% 

% p_Pa/(pc.kB*T_K)
% 
% effective_collision_area = pi*d_m^2
% mfp = 1/(sqrt(2)*effective_collision_area*density)


% Veff = (M_kgPermol/pc.NAvogadro)/rho; % Effective volume
% Vatom = (4/3*pi)/8*Veff;  % Ratio of unit sphere to unit cube
% reff = (Vatom/(4*pi/3))^(1/3)