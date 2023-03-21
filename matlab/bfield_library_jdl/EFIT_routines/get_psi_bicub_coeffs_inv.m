function psi_bicub_coeffs_inv = get_psi_bicub_coeffs_inv(g)
% Calculate bicubic matrix cij such that bicubic interpolation 
%     z(t,u) = sum(i=0:3, sum(j=0:3, cij*ti*uj))
% where
%     t = (x - xi)/(xi+1 - xi)
%     u = (y - yj)/(yj+1 - yj)
% are the relative points on the rectilinear grid. 
% 
% To evaluate z inside of a grid cell, dz/dx, dz/dy, and dz/dxdy are 
% required at each vertex. We use central differences, which will preserve 
% div(B) = 0. The derivatives cannot be calculated for the vertices on the
% boundary, which means that z cannot be evaluated for the boundary cells.
%
% Note: ip_sign = -sign(g.Ip) is ****NOT**** applied to the psirz grid here
%
% Since the bicubic matrix  can be inverted analytically, we use the 
% inverted matrix coefficients directly. 
% 
% JDL

% Inverted matrix coefficients are given in 
% get_bicub_mat_inv()
%      1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
%      0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
%     -3     3     0     0    -2    -1     0     0     0     0     0     0     0     0     0     0
%      2    -2     0     0     1     1     0     0     0     0     0     0     0     0     0     0
%      0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0
%      0     0     0     0     0     0     0     0     0     0     0     0     1     0     0     0
%      0     0     0     0     0     0     0     0    -3     3     0     0    -2    -1     0     0
%      0     0     0     0     0     0     0     0     2    -2     0     0     1     1     0     0
%     -3     0     3     0     0     0     0     0    -2     0    -1     0     0     0     0     0
%      0     0     0     0    -3     0     3     0     0     0     0     0    -2     0    -1     0
%      9    -9    -9     9     6     3    -6    -3     6    -6     3    -3     4     2     2     1
%     -6     6     6    -6    -3    -3     3     3    -4     4    -2     2    -2    -2    -1    -1
%      2     0    -2     0     0     0     0     0     1     0     1     0     0     0     0     0
%      0     0     0     0     2     0    -2     0     0     0     0     0     1     0     1     0
%     -6     6     6    -6    -4    -2     4     2    -3     3    -3     3    -2    -1    -2    -1
%      4    -4    -4     4     2     2    -2    -2     2    -2     2    -2     1     1     1     1
%
% to be applied to 
% [z(0,0), z(1,0), z(0,1), z(1,1), ...
% zx(0,0)dx,  zx(1,0)dx, zx(0,1)dx, zx(1,1)dx, ...
% zy(0,0)dy,  zy(1,0)dy, zy(0,1)dy, zy(1,1)dy, ...
% zxy(0,0)dxdy,  zxy(1,0)dxdy, zxy(0,1)dxdy, zxy(1,1)dxdy]'

%%  Apply sign convention
% psi = g.ip_sign*g.psirz;
psi = g.psirz;

%% Initialize full array, only interior points will be evaluated
nr = g.mw;
nz = g.mh;

ir = 2:nr-1;
iz = 2:nz-1;

dpsidr = nan(nr,nz);
dpsidz = nan(nr,nz);
d2psidrdz = nan(nr,nz);
psi_bicub_coeffs_inv.c00 = nan(nr,nz);
psi_bicub_coeffs_inv.c10 = nan(nr,nz);
psi_bicub_coeffs_inv.c20 = nan(nr,nz);
psi_bicub_coeffs_inv.c30 = nan(nr,nz);
psi_bicub_coeffs_inv.c01 = nan(nr,nz);
psi_bicub_coeffs_inv.c11 = nan(nr,nz);
psi_bicub_coeffs_inv.c21 = nan(nr,nz);
psi_bicub_coeffs_inv.c31 = nan(nr,nz);
psi_bicub_coeffs_inv.c02 = nan(nr,nz);
psi_bicub_coeffs_inv.c12 = nan(nr,nz);
psi_bicub_coeffs_inv.c22 = nan(nr,nz);
psi_bicub_coeffs_inv.c32 = nan(nr,nz);
psi_bicub_coeffs_inv.c03 = nan(nr,nz);
psi_bicub_coeffs_inv.c13 = nan(nr,nz);
psi_bicub_coeffs_inv.c23 = nan(nr,nz);
psi_bicub_coeffs_inv.c33 = nan(nr,nz);

%% Calculate central derivatives on interior points. The 1D derivatives
% would be valid on the boundaries in the orthogonal direction, but we 
% keep them as nan for simplicity.
dpsidr(ir,iz) = (psi(ir+1,iz) - psi(ir-1,iz))/(2*g.dR);
dpsidz(ir,iz) = (psi(ir,iz+1) - psi(ir,iz-1))/(2*g.dZ);
d2psidrdz(ir,iz) = ( psi(ir+1,iz+1) - psi(ir-1,iz+1) ...
                   - psi(ir+1,iz-1) + psi(ir-1,iz-1) )/(4*g.dR*g.dZ);

%% Evaluate coefficients on interior points
psi_bicub_coeffs_inv.c00(ir,iz) = psi(ir,iz);
psi_bicub_coeffs_inv.c10(ir,iz) = dpsidr(ir,iz)*g.dR;
psi_bicub_coeffs_inv.c20(ir,iz) = -3*psi(ir,iz) + 3*psi(ir+1,iz) - 2*dpsidr(ir,iz)*g.dR - dpsidr(ir+1,iz)*g.dR;
psi_bicub_coeffs_inv.c30(ir,iz) =  2*psi(ir,iz) - 2*psi(ir+1,iz) +   dpsidr(ir,iz)*g.dR + dpsidr(ir+1,iz)*g.dR;

psi_bicub_coeffs_inv.c01(ir,iz) = dpsidz(ir,iz)*g.dZ;
psi_bicub_coeffs_inv.c11(ir,iz) = d2psidrdz(ir,iz)*g.dR*g.dZ;
psi_bicub_coeffs_inv.c21(ir,iz) = -3*dpsidz(ir,iz)*g.dZ + 3*dpsidz(ir+1,iz)*g.dZ - 2*d2psidrdz(ir,iz)*g.dR*g.dZ - d2psidrdz(ir+1,iz)*g.dR*g.dZ;
psi_bicub_coeffs_inv.c31(ir,iz) =  2*dpsidz(ir,iz)*g.dZ - 2*dpsidz(ir+1,iz)*g.dZ +   d2psidrdz(ir,iz)*g.dR*g.dZ + d2psidrdz(ir+1,iz)*g.dR*g.dZ;

psi_bicub_coeffs_inv.c02(ir,iz) = -3*psi(ir,iz) + 3*psi(ir,iz+1) - 2*dpsidz(ir,iz)*g.dZ - dpsidz(ir,iz+1)*g.dZ;
psi_bicub_coeffs_inv.c12(ir,iz) = -3*dpsidr(ir,iz)*g.dR + 3*dpsidr(ir,iz+1)*g.dR - 2*d2psidrdz(ir,iz)*g.dR*g.dZ - d2psidrdz(ir,iz+1)*g.dR*g.dZ;
psi_bicub_coeffs_inv.c22(ir,iz) =  9*psi(ir,iz) - 9*psi(ir+1,iz) - 9*psi(ir,iz+1) + 9*psi(ir+1,iz+1) ...
    + 6*dpsidr(ir,iz)*g.dR + 3*dpsidr(ir+1,iz)*g.dR - 6*dpsidr(ir,iz+1)*g.dR - 3*dpsidr(ir+1,iz+1)*g.dR ...
    + 6*dpsidz(ir,iz)*g.dZ - 6*dpsidz(ir+1,iz)*g.dZ + 3*dpsidz(ir,iz+1)*g.dZ - 3*dpsidz(ir+1,iz+1)*g.dZ ...
    + 4*d2psidrdz(ir,iz)*g.dR*g.dZ + 2*d2psidrdz(ir+1,iz)*g.dR*g.dZ + 2*d2psidrdz(ir,iz+1)*g.dR*g.dZ + d2psidrdz(ir+1,iz+1)*g.dR*g.dZ;
psi_bicub_coeffs_inv.c32(ir,iz) = -6*psi(ir,iz) + 6*psi(ir+1,iz) + 6*psi(ir,iz+1) - 6*psi(ir+1,iz+1) ...
    - 3*dpsidr(ir,iz)*g.dR - 3*dpsidr(ir+1,iz)*g.dR + 3*dpsidr(ir,iz+1)*g.dR + 3*dpsidr(ir+1,iz+1)*g.dR ...
    - 4*dpsidz(ir,iz)*g.dZ + 4*dpsidz(ir+1,iz)*g.dZ - 2*dpsidz(ir,iz+1)*g.dZ + 2*dpsidz(ir+1,iz+1)*g.dZ ...
    - 2*d2psidrdz(ir,iz)*g.dR*g.dZ - 2*d2psidrdz(ir+1,iz)*g.dR*g.dZ - d2psidrdz(ir,iz+1)*g.dR*g.dZ - d2psidrdz(ir+1,iz+1)*g.dR*g.dZ;

psi_bicub_coeffs_inv.c03(ir,iz) = 2*psi(ir,iz) - 2*psi(ir,iz+1) + dpsidz(ir,iz)*g.dZ + dpsidz(ir,iz+1)*g.dZ;
psi_bicub_coeffs_inv.c13(ir,iz) = 2*dpsidr(ir,iz)*g.dR - 2*dpsidr(ir,iz+1)*g.dR + d2psidrdz(ir,iz)*g.dR*g.dZ + d2psidrdz(ir,iz+1)*g.dR*g.dZ;
psi_bicub_coeffs_inv.c23(ir,iz) = -6*psi(ir,iz) + 6*psi(ir+1,iz) + 6*psi(ir,iz+1) - 6*psi(ir+1,iz+1) ...
    - 4*dpsidr(ir,iz)*g.dR - 2*dpsidr(ir+1,iz)*g.dR + 4*dpsidr(ir,iz+1)*g.dR + 2*dpsidr(ir+1,iz+1)*g.dR ...
    - 3*dpsidz(ir,iz)*g.dZ + 3*dpsidz(ir+1,iz)*g.dZ - 3*dpsidz(ir,iz+1)*g.dZ + 3*dpsidz(ir+1,iz+1)*g.dZ ...
    - 2*d2psidrdz(ir,iz)*g.dR*g.dZ - d2psidrdz(ir+1,iz)*g.dR*g.dZ - 2*d2psidrdz(ir,iz+1)*g.dR*g.dZ - d2psidrdz(ir+1,iz+1)*g.dR*g.dZ;
psi_bicub_coeffs_inv.c33(ir,iz) =  4*psi(ir,iz) - 4*psi(ir+1,iz) - 4*psi(ir,iz+1) + 4*psi(ir+1,iz+1) ...
    + 2*dpsidr(ir,iz)*g.dR + 2*dpsidr(ir+1,iz)*g.dR - 2*dpsidr(ir,iz+1)*g.dR - 2*dpsidr(ir+1,iz+1)*g.dR ...
    + 2*dpsidz(ir,iz)*g.dZ - 2*dpsidz(ir+1,iz)*g.dZ + 2*dpsidz(ir,iz+1)*g.dZ - 2*dpsidz(ir+1,iz+1)*g.dZ ...
    + d2psidrdz(ir,iz)*g.dR*g.dZ + d2psidrdz(ir+1,iz)*g.dR*g.dZ + d2psidrdz(ir,iz+1)*g.dR*g.dZ + d2psidrdz(ir+1,iz+1)*g.dR*g.dZ;
