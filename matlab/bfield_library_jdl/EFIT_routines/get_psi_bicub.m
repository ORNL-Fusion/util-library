function [psi1,dpsidr,dpsidz] = get_psi_bicub(g,R1,Z1,quiet)
%[psi1,dpsidr,dpsidz] = get_psi_bicub(g,R1,Z1)
% inputs must be row vectors
npts = length(R1);
R1=reshape(R1,1,npts);
Z1=reshape(Z1,1,npts);

if nargin < 4
    quiet = 0;
end

if any(R1 < g.r(1) | R1 > g.r(end) | Z1 < g.z(1) | Z1 > g.z(end))    
%     warning(['Evalution point(s) off grid in get_psi_bicub: ',num2str(R1),num2str(Z1),']'])  
    if ~quiet
        warning('Evalution point(s) off grid in get_psi_bicub')  
    end
    good_inds = find(R1 >= g.r(1) & R1 <= g.r(end) & Z1 >= g.z(1) & Z1 <= g.z(end));
    good_map = (R1 >= g.r(1) & R1 <= g.r(end) & Z1 >= g.z(1) & Z1 <= g.z(end));
else
    good_inds = (1:npts).';
    good_map = good_inds;
end

psi1 = NaN(npts,1);

ir = floor( (R1(good_map) - g.r(1))./g.dR ) + 1;
iz = floor( (Z1(good_map) - g.z(1))./g.dZ ) + 1;

dir(:,1) = (R1(good_map) - g.r(ir))./g.dR;
diz(:,1) = (Z1(good_map) - g.z(iz))./g.dZ;

dir2 = dir.*dir;
dir3 = dir2.*dir;
diz2 = diz.*diz;
diz3 = diz2.*diz;

index = iz + g.mh*(ir-1);
c_bi = g.bicub_coeffs(index,:,:);

psi1(good_inds) = ...
       c_bi(:,1,1)        + c_bi(:,2,1).*dir        + c_bi(:,3,1).*dir2        + c_bi(:,4,1).*dir3       + ...
       c_bi(:,1,2).*diz   + c_bi(:,2,2).*dir.*diz   + c_bi(:,3,2).*dir2.*diz   + c_bi(:,4,2).*dir3.*diz  + ...
       c_bi(:,1,3).*diz2  + c_bi(:,2,3).*dir.*diz2  + c_bi(:,3,3).*dir2.*diz2  + c_bi(:,4,3).*dir3.*diz2 + ...
       c_bi(:,1,4).*diz3  + c_bi(:,2,4).*dir.*diz3  + c_bi(:,3,4).*dir2.*diz3  + c_bi(:,4,4).*dir3.*diz3;

if nargout > 1
%     dspidr = NaN(npts,1);
%     dspidz = NaN(npts,1);
    dpsidr(good_inds) = ...
             (c_bi(:,2,1)        + 2*c_bi(:,3,1).*dir        + 3*c_bi(:,4,1).*dir2        + ...
              c_bi(:,2,2).*diz   + 2*c_bi(:,3,2).*dir.*diz   + 3*c_bi(:,4,2).*dir2.*diz   + ...
              c_bi(:,2,3).*diz2  + 2*c_bi(:,3,3).*dir.*diz2  + 3*c_bi(:,4,3).*dir2.*diz2  + ...
              c_bi(:,2,4).*diz3  + 2*c_bi(:,3,4).*dir.*diz3  + 3*c_bi(:,4,4).*dir2.*diz3)./g.dR ;
    dpsidz(good_inds) = ...
             (  c_bi(:,1,2)        +   c_bi(:,2,2).*dir        +   c_bi(:,3,2).*dir2        +   c_bi(:,4,2).*dir3       + ...
              2*c_bi(:,1,3).*diz   + 2*c_bi(:,2,3).*dir.*diz   + 2*c_bi(:,3,3).*dir2.*diz   + 2*c_bi(:,4,3).*dir3.*diz  + ...
              3*c_bi(:,1,4).*diz2  + 3*c_bi(:,2,4).*dir.*diz2  + 3*c_bi(:,3,4).*dir2.*diz2  + 3*c_bi(:,4,4).*dir3.*diz2)./g.dZ;
end