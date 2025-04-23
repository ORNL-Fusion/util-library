function shape = calc_tokamak_shaping_parameters(R,Z)
% shape = calc_tokamak_shaping_parameters(R,Z)
% R and Z should be arrays of the LCFS shape

% See: T C Luce 2013 Plasma Phys. Control. Fusion 55 095009, 10.1088/0741-3335/55/9/095009

%% Get Rmin and Rmax, Zmin, Zmax perform sanity check
[Rmax,iRmax] = max(R);
[Rmin,iRmin] = min(R);

[Zmax,iZmax] = max(Z);
[Zmin,iZmin] = min(Z);

% Sanity checks
if abs(Z(iRmax)) > (Zmax-Zmin)*0.1
    fprintf('Warning: Rmax used in calc_tokamak_shaping_parameters does not seem to be at the midplane\n')
end
if abs(Z(iRmin)) > (Zmax-Zmin)*0.1
    fprintf('Warning: Rmin used in calc_tokamak_shaping_parameters does not seem to be at the midplane\n')
end
    
%%  minor radius a = 0.5*(max(R) - min(R))
a = 0.5*(Rmax - Rmin);

%% Rgeo
Rgeo = (Rmax + Rmin)/2;

%% Inverse aspect ratio
epsilon = a/Rgeo;

%% Elongation
kappa = (Zmax - Zmin)/(2*a);

%% Triangularity
triUpper = (Rgeo - R(iZmax))/a;
triLower = (Rgeo - R(iZmin))/a;

%% output
shape.Rmin = Rmin;
shape.Rmax = Rmax;
shape.Zmin = Zmin;
shape.Zmax = Zmax;
shape.a = a;
shape.kappa = kappa;
shape.epsilon = epsilon;
shape.triUpper = triUpper;
shape.triLower = triLower;
shape.Rgeo = Rgeo;

end

