function [Br,Bphi,Bz,ierr] = bfield_bs_Bspline(rvec,pvec,zvec,bspline_data,nowarn)

% rr = 0.5;
% zz = -0.2;
% pp = 0.3;
n = length(rvec);
Br = zeros(1,n);
Bphi = zeros(1,n);
Bz = zeros(1,n);
for i=1:n
    rr=rvec(i);
    zz=zvec(i);
    pp=mod(pvec(i),2*pi);
    
    if rr < bspline_data.Rmin || rr > bspline_data.Rmax
        if ~nowarn
            fprintf('Warning!! --> Point off grid in R: R = %f. [Rmin,Rmax] = [%f,%f]\n',rr,bspline_data.Rmin,bspline_data.Rmax);
        end
        ierr = 1; Br = []; Bphi = []; Bz = [];
        return;
    end
    if zz < bspline_data.Zmin || zz > bspline_data.Zmax
        if ~nowarn
            fprintf('Warning!! --> Point off grid in Z: Z = %f. [Zmin,Zmax] = [%f,%f]\n',zz,bspline_data.Zmin,bspline_data.Zmax);
        end
        ierr = 1; Br = []; Bphi = []; Bz = [];
        return;
    end    
    

    dAr_dz = dbs3dr(0,1,0,rr,zz,pp,bspline_data.nord,bspline_data.nord,bspline_data.nord,bspline_data.Rnot,bspline_data.Znot,bspline_data.Phinot,bspline_data.nr,bspline_data.nz,bspline_data.nphi,bspline_data.Arcoeff);
    dAr_dp = dbs3dr(0,0,1,rr,zz,pp,bspline_data.nord,bspline_data.nord,bspline_data.nord,bspline_data.Rnot,bspline_data.Znot,bspline_data.Phinot,bspline_data.nr,bspline_data.nz,bspline_data.nphi,bspline_data.Arcoeff);
    dAz_dr = dbs3dr(1,0,0,rr,zz,pp,bspline_data.nord,bspline_data.nord,bspline_data.nord,bspline_data.Rnot,bspline_data.Znot,bspline_data.Phinot,bspline_data.nr,bspline_data.nz,bspline_data.nphi,bspline_data.Azcoeff);
    dAz_dp = dbs3dr(0,0,1,rr,zz,pp,bspline_data.nord,bspline_data.nord,bspline_data.nord,bspline_data.Rnot,bspline_data.Znot,bspline_data.Phinot,bspline_data.nr,bspline_data.nz,bspline_data.nphi,bspline_data.Azcoeff);
    dAp_dr = dbs3dr(1,0,0,rr,zz,pp,bspline_data.nord,bspline_data.nord,bspline_data.nord,bspline_data.Rnot,bspline_data.Znot,bspline_data.Phinot,bspline_data.nr,bspline_data.nz,bspline_data.nphi,bspline_data.Aphicoeff);
    dAp_dz = dbs3dr(0,1,0,rr,zz,pp,bspline_data.nord,bspline_data.nord,bspline_data.nord,bspline_data.Rnot,bspline_data.Znot,bspline_data.Phinot,bspline_data.nr,bspline_data.nz,bspline_data.nphi,bspline_data.Aphicoeff);
    
    % Ar   = dbs3dr(0,0,0,rr,zz,pp,bspline_data.nord,bspline_data.nord,bspline_data.nord,bspline_data.Rnot,bspline_data.Znot,bspline_data.Phinot,bspline_data.nr,bspline_data.nz,bspline_data.nphi,bspline_data.Arcoeff);
    % Az   = dbs3dr(0,0,0,rr,zz,pp,bspline_data.nord,bspline_data.nord,bspline_data.nord,bspline_data.Rnot,bspline_data.Znot,bspline_data.Phinot,bspline_data.nr,bspline_data.nz,bspline_data.nphi,bspline_data.Azcoeff);
    Aphi = dbs3dr(0,0,0,rr,zz,pp,bspline_data.nord,bspline_data.nord,bspline_data.nord,bspline_data.Rnot,bspline_data.Znot,bspline_data.Phinot,bspline_data.nr,bspline_data.nz,bspline_data.nphi,bspline_data.Aphicoeff);
    
    Br(i) = 1/rr * dAz_dp - dAp_dz;
    Bphi(i) = dAr_dz - dAz_dr;
    Bz(i) = 1/rr *Aphi + dAp_dr - 1/rr*dAr_dp;
    ierr = 0;
end