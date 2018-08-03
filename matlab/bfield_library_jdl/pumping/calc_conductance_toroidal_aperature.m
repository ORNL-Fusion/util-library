clearvars;
% See Livesey

% Define aperature (Here assumed vertical opening at constant Rmajor)
H_ap = 0.035; % Aperature height (m)
R_ap = 1.65;  % Major radius of aperature (m)

% Define duct (also assumed to be in Rmajor direction)
L_duct = 0.1; % Duct length (m)
nL = 100; % Number of points for integration along duct

% Define gas parameters
T = 300;  % Temperature (K)
A = 2.014; % Atomic mass of constituent atom (u)
Nmolec = 2; % Number of atoms per molecule

% Constants
amu = 1.66e-27;   % kg
R0  = 8.3144598;  % Universal gas constant (J/mol/K)
Na  = 6.022e23;   % Avogadro's number 
kB  = 1.38e-23;   % Boltzmann constant (J/K)


% Calculate
Mm = Nmolec*A*amu*Na;  % Molar mass  (kg/mol)
Q1D = sqrt(R0*T/(2*pi*Mm)); % 1D flow rate (m/s)
Area = 2*pi*R_ap*H_ap;  
Factor = 2*pi*Q1D  %  C = Factor*R_ap*H_ap  ~ 1972 for D2
C_ap = Factor *R_ap*H_ap  % Conductance (m^3/s)

M_molec = amu*Nmolec*A;
v_mean = sqrt(8*kB*T/M_molec/pi);
Q1D_v2 = v_mean/4;
Factor_v2 = 2*pi*Q1D_v2

% Now consider a duct
R_duct = R_ap + linspace(0,L_duct,nL);
dR = R_duct(2) - R_duct(1);

% Knudsen formula for long duct (cylindrical)
perim_duct = 2*2*pi*R_duct;
area_duct  = 2*pi*H_ap*R_duct; 
factor_duct = 1/(sum(perim_duct./area_duct.^2)*dR);
C_m_Knud = 4/3*v_mean * factor_duct
factor_duct_v2 = H_ap^2*pi/log((R_ap+L_duct)/R_ap)  % Analytic integral of above

% Smoluchowski formula
% f
% C_m_Smol = v

% Total conductance (series formula)
Ctot = 1 / ( (1/C_m_Knud) + (1/C_ap) )

% particle density per Pascal
P_Pa = 1;
density_per_Pa = P_Pa*Na/(R0*T)

