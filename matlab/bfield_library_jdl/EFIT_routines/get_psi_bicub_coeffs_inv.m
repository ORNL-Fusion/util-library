function psi_bicub_coeffs = get_psi_bicub_coeffs_inv(g)

nr = g.mw;
nz = g.mh;

psi = g.ip_sign*g.psirz;

% Calculate derivatives
dpsidr = (circshift(psi,[-1,0])-circshift(psi,[1,0]))/(2*g.dR);
dpsidz = (circshift(psi,[0,-1])-circshift(psi,[0,1]))/(2*g.dZ);
d2psidrdz = (circshift(psi,[-1,-1])-circshift(psi,[1,-1])-circshift(psi,[-1,1])+circshift(psi,[1,1]))/(4*g.dR*g.dZ);


ir = 1:nr-1;
iz = 1:nz-1;

% fix this indexing
psi_bicub_coeffs.c00 = psi(ir,iz);
psi_bicub_coeffs.c10 = dpsidr(ir,iz)*g.dR;
psi_bicub_coeffs.c20 = -3*psi(ir,iz) + 3*psi(ir+1,iz) - 2*dpsidr(ir,iz)*g.dR - dpsidr(ir+1,iz)*g.dR;
psi_bicub_coeffs.c30 =  2*psi(ir,iz) - 2*psi(ir+1,iz) +   dpsidr(ir,iz)*g.dR + dpsidr(ir+1,iz)*g.dR;

psi_bicub_coeffs.c01 = dpsidz(ir,iz)*g.dZ;
psi_bicub_coeffs.c11 = d2psidrdz(ir,iz)*g.dR*g.dZ;
psi_bicub_coeffs.c21 = -3*dpsidz(ir,iz)*g.dZ + 3*dpsidz(ir+1,iz)*g.dZ - 2*d2psidrdz(ir,iz)*g.dR*g.dZ - d2psidrdz(ir+1,iz)*g.dR*g.dZ;
psi_bicub_coeffs.c31 =  2*dpsidz(ir,iz)*g.dZ - 2*dpsidz(ir+1,iz)*g.dZ +   d2psidrdz(ir,iz)*g.dR*g.dZ + d2psidrdz(ir+1,iz)*g.dR*g.dZ;

psi_bicub_coeffs.c02 = -3*psi(ir,iz) + 3*psi(ir,iz+1) - 2*dpsidz(ir,iz)*g.dZ - dpsidz(ir,iz+1)*g.dZ;
psi_bicub_coeffs.c12 = -3*dpsidr(ir,iz)*g.dR + 3*dpsidr(ir,iz+1)*g.dR - 2*d2psidrdz(ir,iz)*g.dR*g.dZ - d2psidrdz(ir,iz+1)*g.dR*g.dZ;
psi_bicub_coeffs.c22 =  9*psi(ir,iz) - 9*psi(ir+1,iz) - 9*psi(ir,iz+1) + 9*psi(ir+1,iz+1) ...
    + 6*dpsidr(ir,iz)*g.dR + 3*dpsidr(ir+1,iz)*g.dR - 6*dpsidr(ir,iz+1)*g.dR - 3*dpsidr(ir+1,iz+1)*g.dR ...
    + 6*dpsidz(ir,iz)*g.dZ - 6*dpsidz(ir+1,iz)*g.dZ + 3*dpsidz(ir,iz+1)*g.dZ - 3*dpsidz(ir+1,iz+1)*g.dZ ...
    + 4*d2psidrdz(ir,iz)*g.dR*g.dZ + 2*d2psidrdz(ir+1,iz)*g.dR*g.dZ + 2*d2psidrdz(ir,iz+1)*g.dR*g.dZ + d2psidrdz(ir+1,iz+1)*g.dR*g.dZ;
psi_bicub_coeffs.c32 = -6*psi(ir,iz) + 6*psi(ir+1,iz) + 6*psi(ir,iz+1) - 6*psi(ir+1,iz+1) ...
    - 3*dpsidr(ir,iz)*g.dR - 3*dpsidr(ir+1,iz)*g.dR + 3*dpsidr(ir,iz+1)*g.dR + 3*dpsidr(ir+1,iz+1)*g.dR ...
    - 4*dpsidz(ir,iz)*g.dZ + 4*dpsidz(ir+1,iz)*g.dZ - 2*dpsidz(ir,iz+1)*g.dZ + 2*dpsidz(ir+1,iz+1)*g.dZ ...
    - 2*d2psidrdz(ir,iz)*g.dR*g.dZ - 2*d2psidrdz(ir+1,iz)*g.dR*g.dZ - d2psidrdz(ir,iz+1)*g.dR*g.dZ - d2psidrdz(ir+1,iz+1)*g.dR*g.dZ;

psi_bicub_coeffs.c03 = 2*psi(ir,iz) - 2*psi(ir,iz+1) + dpsidz(ir,iz)*g.dZ + dpsidz(ir,iz+1)*g.dZ;
psi_bicub_coeffs.c13 = 2*dpsidr(ir,iz)*g.dR - 2*dpsidr(ir,iz+1)*g.dR + d2psidrdz(ir,iz)*g.dR*g.dZ + d2psidrdz(ir,iz+1)*g.dR*g.dZ;
psi_bicub_coeffs.c23 = -6*psi(ir,iz) + 6*psi(ir+1,iz) + 6*psi(ir,iz+1) - 6*psi(ir+1,iz+1) ...
    - 4*dpsidr(ir,iz)*g.dR - 2*dpsidr(ir+1,iz)*g.dR + 4*dpsidr(ir,iz+1)*g.dR + 2*dpsidr(ir+1,iz+1)*g.dR ...
    - 3*dpsidz(ir,iz)*g.dZ + 3*dpsidz(ir+1,iz)*g.dZ - 3*dpsidz(ir,iz+1)*g.dZ + 3*dpsidz(ir+1,iz+1)*g.dZ ...
    - 2*d2psidrdz(ir,iz)*g.dR*g.dZ - d2psidrdz(ir+1,iz)*g.dR*g.dZ - 2*d2psidrdz(ir,iz+1)*g.dR*g.dZ - d2psidrdz(ir+1,iz+1)*g.dR*g.dZ;
psi_bicub_coeffs.c33 =  4*psi(ir,iz) - 4*psi(ir+1,iz) - 4*psi(ir,iz+1) + 4*psi(ir+1,iz+1) ...
    + 2*dpsidr(ir,iz)*g.dR + 2*dpsidr(ir+1,iz)*g.dR - 2*dpsidr(ir,iz+1)*g.dR - 2*dpsidr(ir+1,iz+1)*g.dR ...
    + 2*dpsidz(ir,iz)*g.dZ - 2*dpsidz(ir+1,iz)*g.dZ + 2*dpsidz(ir,iz+1)*g.dZ - 2*dpsidz(ir+1,iz+1)*g.dZ ...
    + d2psidrdz(ir,iz)*g.dR*g.dZ + d2psidrdz(ir+1,iz)*g.dR*g.dZ + d2psidrdz(ir,iz+1)*g.dR*g.dZ + d2psidrdz(ir+1,iz+1)*g.dR*g.dZ;
