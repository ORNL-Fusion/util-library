function [Br,Bz,Bphi,ierr]=bfield_Aspline(R,Z,phi_radian,Arcoeff,Azcoeff,Aphicoeff,spline_info,nowarn)
if nargin < 8
    nowarn = 0;
end

npts = length(R);
Br=zeros(npts,1);
Bz=zeros(npts,1);
Bphi=zeros(npts,1);

nord = spline_info.nord;
Rnot = spline_info.Rnot;
Znot = spline_info.Znot;
Phinot = spline_info.Phinot;
nr = spline_info.nr;
nz = spline_info.nz;
nphi = spline_info.nphi;
nc = length(spline_info.scale_factor);

ierr = 0;
for i = 1:npts
    rr = R(i);
    zz = Z(i);
    pp = phi_radian(i);
    if ~nowarn && (rr < spline_info.Rmin || rr > spline_info.Rmax)
        warning(['Point off grid in R: R = ',num2str(rr),'. [Rmin,Rmax] = [',num2str(spline_info.Rmin),',',num2str(spline_info.Rmax),']'])
        ierr = 1;
    end
    if ~nowarn && (zz < spline_info.Zmin || zz > spline_info.Zmax)
        warning(['Point off grid in Z: Z = ',num2str(zz),'. [Zmin,Zmax] = [',num2str(spline_info.Zmin),',',num2str(spline_info.Zmax),']'])
        ierr = 1;   
    end
    
    pp = mod(pp,2*pi/spline_info.nsym);
    % Br = 0;
    % Bphi = 0;
    % Bz = 0;
    

    
    Br_tmp = zeros(nc,1);
    Bphi_tmp = zeros(nc,1);
    Bz_tmp = zeros(nc,1);
    for ic = 1:nc
        Br_tmp = zeros(nc,1);
        Bphi_tmp = zeros(nc,1);
        Bz_tmp = zeros(nc,1);
        dAr_dz = dbs3dr(0,1,0,rr,zz,pp,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Arcoeff{ic});
        dAr_dp = dbs3dr(0,0,1,rr,zz,pp,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Arcoeff{ic});
        dAz_dr = dbs3dr(1,0,0,rr,zz,pp,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Azcoeff{ic});
        dAz_dp = dbs3dr(0,0,1,rr,zz,pp,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Azcoeff{ic});
        dAp_dr = dbs3dr(1,0,0,rr,zz,pp,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Aphicoeff{ic});
        dAp_dz = dbs3dr(0,1,0,rr,zz,pp,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Aphicoeff{ic});
        
        Aphi = dbs3vl(rr,zz,pp,nord,nord,nord,Rnot,Znot,Phinot,nr,nz,nphi,Aphicoeff{ic});
        
        Br_tmp(ic) = 1/rr * dAz_dp - dAp_dz;
        Bphi_tmp(ic) = dAr_dz - dAz_dr;
        Bz_tmp(ic) = 1/rr *Aphi + dAp_dr - 1/rr*dAr_dp;
    end
    
    Br(i) = sum(Br_tmp.*spline_info.scale_factor);
    Bz(i) = sum(Bz_tmp.*spline_info.scale_factor);
    Bphi(i) = sum(Bphi_tmp.*spline_info.scale_factor);
    
end