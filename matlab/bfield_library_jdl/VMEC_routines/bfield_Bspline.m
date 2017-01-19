function [Br,Bz,Bphi]=bfield_Bspline(rr,zz,pp,Brcoeff,Bzcoeff,Bphicoeff,spline_info)   

pp = mod(pp,2*pi/spline_info.nsym);
% Br = 0;
% Bphi = 0;
% Bz = 0;

nord = spline_info.nord;
Rnot = spline_info.Rnot;
Znot = spline_info.Znot;
Phinot = spline_info.Phinot;
nr = spline_info.nr;
nz = spline_info.nz;
nphi = spline_info.nphi;
nc = length(spline_info.scale_factor);

for ic = 1:nc
    Br(ic)   = dbs3vl(rr,zz,pp,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Brcoeff{ic});
    Bz(ic)   = dbs3vl(rr,zz,pp,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Bzcoeff{ic});
    Bphi(ic) = dbs3vl(rr,zz,pp,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Bphicoeff{ic});       
end

Br = sum(Br.*spline_info.scale_factor);
Bz = sum(Bz.*spline_info.scale_factor);
Bphi = sum(Bphi.*spline_info.scale_factor);

