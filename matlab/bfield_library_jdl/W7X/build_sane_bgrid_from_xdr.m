% function build_sane_bgrid_from_xdr(b)

out_path = 'C:\Work\Stellarator\ALL_W7X_WORK\xdr_dump_read\OUTPUT\';
fname = 'field181x181x96.w7x.1000_1000_1000_1000_+0750_+0750.vac.out';
xdr = read_xdr_dump_file(out_path,fname);

rmin = xdr.rnull - xdr.ronull;
dr = xdr.ronull/xdr.knull;
delta_R = (xdr.k2-1)*dr;
rmax = rmin + delta_R;

dz = xdr.eta*dr;
zmin = -(xdr.k2-1)*dz/2;  % assume zonull = 0
zmax = -zmin;

disp(['Grid RZ domain [Rmin,Rmax,Zmin,Zmax] = ',num2str([rmin,rmax,zmin,zmax])])

% info for full field period
phimin_fp = 0;
phimax_fp = 2*pi/xdr.nperio;
delta_phi_fp = phimax_fp - phimin_fp;
dphi_fp = delta_phi_fp/xdr.ialfa;

% for half field period
phimin_hfp = 0;
phimax_hfp = pi/xdr.nperio;
delta_phi_hfp = phimax_hfp - phimin_hfp;
dphi_hfp = delta_phi_hfp/(xdr.iald21 - 1);

% build hfp grid
nphi_hfp = xdr.iald21;
nphi_fp = xdr.ialfa + 1;
nr = xdr.k2;
nz = xdr.k2;

phi_array = linspace(phimin_hfp,phimax_hfp,nphi_hfp);
% phi_array = linspace(phimin_fp,phimax_fp,nphi_fp);
r_array = linspace(rmin,rmax,nr);
z_array = linspace(zmin,zmax,nz);

[bgrid.rg,bgrid.pg,bgrid.zg] = meshgrid(r_array,phi_array,z_array);

bgrid.brg = xdr.brg;
bgrid.bfg = xdr.bfg;
bgrid.bzg = xdr.bzg;

bgrid.stellsym = 1;
bgrid.nfp = xdr.nperio;

bgrid.rmin = rmin;
bgrid.rmax = rmax;
bgrid.zmin = zmin;
bgrid.zmax = zmax;


save([fname(1:end-3),'bgrid.mat'],'bgrid')


