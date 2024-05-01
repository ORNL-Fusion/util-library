function pc = phys_const(want_name)
% If a name is provided, only that variable is returned

% M_m         = 4.028;
% M_a         = 2.014;

% pc.amu  = 1.660539040e-27;     % kg
pc.me   = 9.1093837015e-31;      % kg
pc.mp   = 1.67262192369e-27;     % kg
pc.e    = 1.602176634e-19;     % C
pc.eV   = 1.602176634e-19;     % C
pc.eps0 = 8.854187817e-12;     % F/m
pc.mu0  = 4*pi*1e-7;           % N/A**2
pc.c    = 2.99792458e8;        % m/s
pc.kB   = 1.38064852e-23;      % J/K
pc.convert.Pa_to_mTorr = 7.5;
pc.convert.K_to_eV     = pc.kB/pc.eV;
pc.convert.eV_to_K     = pc.eV/pc.kB;
pc.convert.C_to_K_Add  = 273.15;
pc.convert.L_to_m3 = 1e-3;
pc.convert.m3_to_L = 1e3;
pc.NAvogadro = 6.02214076e+23; % 1/mol
pc.Rgas = pc.NAvogadro*pc.kB;
pc.amu0 = 1/pc.NAvogadro/1000; % kg

if nargin == 1
    try
        pc = pc.(want_name);
    catch
        error('Did not recognize name %s',want_name)
    end
end