function psi_bicub_coeffs_inv = get_psi_bicub_coeffs_inv(g)
% Calculate bicubic matrix cij such that bicubic interpolation
%     z(t,u) = sum(i=0:3, sum(j=0:3, cij*ti*uj))
% where
%     t = (x - xi)/(xi+1 - xi)
%     u = (y - yj)/(yj+1 - yj)
% are the relative points on the rectilinear grid.
%
% To evaluate z inside of a grid cell, dz/dx, dz/dy, and dz/dxdy are
% required at each vertex. We use central differences internally and
% second-order one-sided differences on the grid edges, allowing all cells to
% be evaluated while preserving div(B) = 0.
%
% Note: ip_sign = -sign(g.Ip) is ****NOT**** applied to the psirz grid here
%
% Since the bicubic matrix  can be inverted analytically, we use the
% inverted matrix coefficients directly.
%
% JDL

ADD_EDGES = 1; % Use second-order one-sided differences on grid edges.

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

%% Initialize derivative and coefficient arrays
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

if ADD_EDGES
    dpsidr(2:nr-1,:) = (psi(3:nr,:) - psi(1:nr-2,:))/(2*g.dR);
    dpsidr(1,:) = (-3*psi(1,:) + 4*psi(2,:) - psi(3,:))/(2*g.dR);
    dpsidr(nr,:) = (3*psi(nr,:) - 4*psi(nr-1,:) + psi(nr-2,:))/(2*g.dR);

    dpsidz(:,2:nz-1) = (psi(:,3:nz) - psi(:,1:nz-2))/(2*g.dZ);
    dpsidz(:,1) = (-3*psi(:,1) + 4*psi(:,2) - psi(:,3))/(2*g.dZ);
    dpsidz(:,nz) = (3*psi(:,nz) - 4*psi(:,nz-1) + psi(:,nz-2))/(2*g.dZ);

    d2psidrdz(:,2:nz-1) = (dpsidr(:,3:nz) - dpsidr(:,1:nz-2))/(2*g.dZ);
    d2psidrdz(:,1) = (-3*dpsidr(:,1) + 4*dpsidr(:,2) - dpsidr(:,3))/(2*g.dZ);
    d2psidrdz(:,nz) = (3*dpsidr(:,nz) - 4*dpsidr(:,nz-1) + dpsidr(:,nz-2))/(2*g.dZ);

    ir = 1:nr-1;
    iz = 1:nz-1;
else
    dpsidr(ir,iz) = (psi(ir+1,iz) - psi(ir-1,iz))/(2*g.dR);
    dpsidz(ir,iz) = (psi(ir,iz+1) - psi(ir,iz-1))/(2*g.dZ);
    d2psidrdz(ir,iz) = (psi(ir+1,iz+1) - psi(ir-1,iz+1) ...
        - psi(ir+1,iz-1) + psi(ir-1,iz-1))/(4*g.dR*g.dZ);
end

%% Evaluate coefficients on valid cells
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
