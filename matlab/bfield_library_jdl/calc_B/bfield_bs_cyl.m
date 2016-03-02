function [Br,Bphi,Bz,Btot]=bfield_bs_cyl(P_r,P_phi,P_z,coil,current,nowarn)
% [Br,Bphi,Bz]=bfield_bs_cyl(P_r,P_phi,P_z,coil,current,nowarn)
if nargin < 6
    nowarn = 0;
end
    
if isrow(P_r)
    P_r = P_r.';
    P_phi = P_phi.';
    P_z = P_z.';
end

cp = cos(P_phi);
sp = sin(P_phi);

P_x = P_r.*cp;
P_y = P_r.*sp;
[Bx,By,Bz]=bfield_bs_jdl(P_x,P_y,P_z,coil,current);

Br = Bx.*cp + By.*sp;
Bphi = -Bx.*sp + By.*cp;

if nargout > 3
    Btot = sqrt(Br.^2 + Bphi.^2 + Bz.^2);
end