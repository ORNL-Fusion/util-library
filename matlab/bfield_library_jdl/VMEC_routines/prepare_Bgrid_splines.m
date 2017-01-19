function [Brcoeff,Bzcoeff,Bphicoeff,spline_info] = prepare_Bgrid_splines(mgrid)

nord = 5;


nr = mgrid.nr;
nz = mgrid.nz;
nphi = mgrid.nphi + 1;

R =  mgrid.R;
Z =  mgrid.Z;
phi= mgrid.phi; phi(nphi) = 2*pi/mgrid.nsym;

% Prepare grid nodes
Rnot   = dbsnak(nr,R,nord);
Znot   = dbsnak(nz,Z,nord);
Phinot = dbsnak(nphi,phi,nord);

tic
fprintf('Preparing B-spline coefficients for B\n')
tic;
nc = length(mgrid.scale_factor);
for ic = 1:nc        
    fprintf('Working on ic %d of %d\n',ic,nc)
    Br   = mgrid.Br{ic};
    Bz   = mgrid.Bz{ic};
    Bphi = mgrid.Bphi{ic};
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
spline_info.scale_factor = mgrid.scale_factor;
spline_info.nsym = mgrid.nsym;