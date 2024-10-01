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
    g.aminor = (max(g.bdry(1,:))-min(g.bdry(1,:)))/2;    
end

g.nGW = [];
if ~isempty(g.cpasma) && ~isempty(g.aminor) 
    g.nGW = g.cpasma/1e6/pi/g.aminor^2*1e20;
end

if isempty(g.cpasma)
    g.ip_sign = 1;
else
    g.ip_sign = -sign(g.cpasma);
end

%% Flux interpolation coefficients
g.bicub_coeffs_inv = get_psi_bicub_coeffs_inv(g);

%% Fit F
g.fpol_coeffs = polyfit(g.pn,g.fpol,7);

end