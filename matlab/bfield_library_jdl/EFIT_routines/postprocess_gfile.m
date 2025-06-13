function g = postprocess_gfile(g)
%
% Postprocessing
%
g.dR = g.xdim/(g.mw-1);
g.dZ = g.zdim/(g.mh-1);

g.r = zeros(1,g.mw);
g.z = zeros(1,g.mh);
g.pn = zeros(1,g.mw);
for i = 0:g.mw-1
    g.r(1,i+1) = g.rgrid1 + g.dR*i;
    g.pn(1,i+1) = i/(g.mw-1);
end

for i = 0:g.mh-1
    g.z(1,i+1) = g.zmid - 0.5*g.zdim + g.dZ*i;
end

g.aminor = [];
if ~isempty(g.bdry)
    shape = calc_tokamak_shaping_parameters(g.bdry(1,:),g.bdry(2,:));
    g.Rmin = shape.Rmin;
    g.Rmax = shape.Rmax;
    g.Zmin = shape.Zmin;
    g.Zmax = shape.Zmax;
    g.aminor = shape.a;
    g.kappa = shape.kappa;
    g.epsilon = shape.epsilon;
    g.triUpper = shape.triUpper;
    g.triLower = shape.triLower;
    g.Rgeo = shape.Rgeo;
end

g.nGW = [];
if ~isempty(g.cpasma) && ~isempty(g.aminor) 
    g.nGW = g.cpasma/1e6/pi/g.aminor^2*1e20;
end

if isempty(g.cpasma)
    g.ip_sign = 1;
else
    if abs(g.cpasma) < 1e-6
        g.ip_sign = 1;
    else
        g.ip_sign = -sign(g.cpasma);
    end
end

%% Flux interpolation coefficients
g.bicub_coeffs_inv = get_psi_bicub_coeffs_inv(g);

%% Fit F
g.fpol_coeffs = polyfit(g.pn,g.fpol,7);

end