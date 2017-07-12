function [Br,Bz,Bphi,ierr] = bfield_xdr(R,Z,P_rad,xdr,nowarn)

npts = length(R);
Br=zeros(npts,1);
Bz=zeros(npts,1);
Bphi=zeros(npts,1);

ierr_tmp = zeros(npts,1);

for i = 1:npts

    [Br(i,1),Bz(i,1),Bphi(i,1),ierr_tmp(i)] = bint_xdr(R(i),Z(i),P_rad(i),xdr,nowarn);


end

ierr = any(ierr_tmp);

