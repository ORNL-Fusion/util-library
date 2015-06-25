function [rs,phis,zs] = symmetrize_vec(r,phi,z,nfp,stellsym)
%this function takes an input position in r,phi,z space
%which can be anywhere in the stellarator and returns
%a position poss which is equivalent but lies in the
%first period (in the case of no "stellarator symmetry"
%or in the first half field period if there is 
%"stellarator symmetry"
phi = mod(phi,2*pi/nfp);
if stellsym
    i = find(phi > pi/nfp);
    z(i) = -z(i);
    phi(i) = 2*pi/nfp - phi(i);
end

zs = z;
rs = r;
phis = phi;