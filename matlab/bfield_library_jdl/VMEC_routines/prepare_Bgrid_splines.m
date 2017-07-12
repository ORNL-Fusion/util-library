function [Brcoeff,Bzcoeff,Bphicoeff,spline_info] = prepare_Bgrid_splines(mgrid_or_bmw)

nord = 5;

% Check that last plane is missing
if abs(diff(mod(mgrid_or_bmw.phi([1,end]),2*pi/mgrid_or_bmw.nsym))) < 2*eps
    error('This routine assumes that the grid does not include the last periodic plane, would be easy to generalize')
end

nr = mgrid_or_bmw.nr;
nz = mgrid_or_bmw.nz;
nphi = mgrid_or_bmw.nphi + 1;

R =  mgrid_or_bmw.R;
Z =  mgrid_or_bmw.Z;
phi= mgrid_or_bmw.phi; phi(nphi) = 2*pi/mgrid_or_bmw.nsym;

% Prepare grid nodes
Rnot   = dbsnak(nr,R,nord);
Znot   = dbsnak(nz,Z,nord);
Phinot = dbsnak(nphi,phi,nord);

tic
fprintf('Preparing B-spline coefficients for B\n')
tic;
nc = length(mgrid_or_bmw.scale_factor);
for ic = 1:nc        
    fprintf('Working on ic %d of %d\n',ic,nc)
    Br   = mgrid_or_bmw.Br{ic};
    Bz   = mgrid_or_bmw.Bz{ic};
    Bphi = mgrid_or_bmw.Bphi{ic};
    Br(:,:,nphi)   = Br(:,:,1);    
    Bz(:,:,nphi)   = Bz(:,:,1);    
    Bphi(:,:,nphi) = Bphi(:,:,1);

    Brcoeff{ic}  = dbs3in(nr,R,nz,Z,nphi,phi,Br  ,nr,nz,nord,nord,nord,Rnot,Znot,Phinot);
    Bzcoeff{ic}  = dbs3in(nr,R,nz,Z,nphi,phi,Bz  ,nr,nz,nord,nord,nord,Rnot,Znot,Phinot);
    Bphicoeff{ic}= dbs3in(nr,R,nz,Z,nphi,phi,Bphi,nr,nz,nord,nord,nord,Rnot,Znot,Phinot);
end
fprintf('Spline prep took %f seconds.\n',toc); tic;

spline_info.nord = nord;
spline_info.R = R;
spline_info.Z = Z;
spline_info.phi = phi;
spline_info.nr = nr;
spline_info.nz = nz;
spline_info.nphi = nphi;
spline_info.Rnot = Rnot;
spline_info.Znot = Znot;
spline_info.Phinot = Phinot;
spline_info.scale_factor = mgrid_or_bmw.scale_factor;
spline_info.nsym = mgrid_or_bmw.nsym;
spline_info.Rmin = mgrid_or_bmw.rmin;
spline_info.Rmax = mgrid_or_bmw.rmax;
spline_info.Zmin = mgrid_or_bmw.zmin;
spline_info.Zmax = mgrid_or_bmw.zmax;