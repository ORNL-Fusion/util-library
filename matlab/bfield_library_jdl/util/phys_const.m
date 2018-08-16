function pc = phys_const

% M_m         = 4.028;
% M_a         = 2.014;

pc.amu  = 1.660539040e-27;     % kg
pc.me   = 9.10938215e-31;      % kg
pc.mp   = 1.672621637e-27;     % kg
pc.e    = 1.602176487e-19;     % C
pc.eV   = 1.602176487e-19;     % C
pc.eps0 = 8.854187817e-12;     % F/m
pc.mu0  = 4*pi*1e-7;           % N/A**2
pc.c    = 2.99792458e8;        % m/s
pc.kB   = 1.38064852e-23;      % J/K
pc.convert.Pa_to_mTorr = 7.5;
pc.convert.K_to_eV     = pc.kB/pc.eV;
pc.convert.eV_to_K     = pc.eV/pc.kB;