function [dpsidr,dpsidz,ierr] = calc_psi_derivs(g,R,Z,nowarn)
% -----------------------------------------------------------------------
% Get psi derivatives at a point from a gfile
% psi = -sign(Ip)*psirz(R,Z), where psirz is the value in the gfile
%
% Uses bicubic interpoation based on psirz grid. 
% -----------------------------------------------------------------------
if nargin < 4
    nowarn = 0;
end

npts = length(R);
R=reshape(R,1,npts);
Z=reshape(Z,1,npts);

[ir,iz,index,ierr] = calc_bicub_interpolation_inds(g,R,Z,nowarn);

% Handle NaN points 
% Assign index of 1 for now -- they will get marked as nan below
indNaN = isnan(ir);
ir(indNaN) = 1;
iz(indNaN) = 1;
index(indNaN) = 1;

dr = (R - g.r(ir))/g.dR;
dz = (Z - g.z(iz))/g.dZ;

c = g.bicub_coeffs_inv;

c10 = c.c10(index);
c20 = c.c20(index);
c30 = c.c30(index);

c01 = c.c01(index);
c11 = c.c11(index);
c21 = c.c21(index);
c31 = c.c31(index);

c02 = c.c02(index);
c12 = c.c12(index);
c22 = c.c22(index);
c32 = c.c32(index);

c03 = c.c03(index);
c13 = c.c13(index);
c23 = c.c23(index);
c33 = c.c33(index);

drr = dr.*dr;
drrr = drr.*dr;
dzz = dz.*dz;
dzzz = dzz.*dz;

drz   = dr.*dz;
drzz  = dr.*dzz;
drzzz = dr.*dzzz;
drrz  = drr.*dz;
drrzz = drr.*dzz;
drrzzz = drr.*dzzz;
drrrz  = drrr.*dz;
drrrzz = drrr.*dzz;


dpsidr = g.ip_sign*(...
    c10       + 2*c20.*dr    + 3*c30.*drr   + ...
    c11.*dz   + 2*c21.*drz   + 3*c31.*drrz  + ...
    c12.*dzz  + 2*c22.*drzz  + 3*c32.*drrzz + ...
    c13.*dzzz + 2*c23.*drzzz + 3*c33.*drrzzz)./g.dR;

dpsidz = g.ip_sign*(...
    c01      +   c11.*dr   +   c21.*drr   +   c31.*drrr   + ...
    2*c02.*dz  + 2*c12.*drz  + 2*c22.*drrz  + 2*c32.*drrrz  + ...
    3*c03.*dzz + 3*c13.*drzz + 3*c23.*drrzz + 3*c33.*drrrzz)./g.dZ;

dpsidr(ierr) = nan;
dpsidz(ierr) = nan;
