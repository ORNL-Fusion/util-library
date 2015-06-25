function [bval,idiv] = bint_xdr(b,xvec)
% clearvars;
% xvec = [6.0,0,0];
% out_path = 'C:\Work\Stellarator\ALL_W7X_WORK\xdr_dump_read\OUTPUT\';
% fname = 'field181x181x96.w7x.1000_1000_1000_1000_+0750_+0750.vac.out';
% b = read_xdr_dump_file(out_path,fname);


%  input : xvec(3), where
% 
%       xvec(1)  = R          xvec(2) = phi         xvec(3) = z
% 
%  output: bval(3), where
% 
%       bval(1)  = br         bval(2)  = bphi       bval(3)  = bz
%       idiv =   0     if bmod.ne.0 and xvec inside of int.region
%                1     if bmod.eq.0
%                   or if xvec is out of the interpolation region
% ***********************************************************************

idebug = 0;
if idebug
    rmin = b.rnull - b.ronull;
    dr = b.ronull/b.knull;
    delta_R = (b.k2-1)*dr;
    rmax = rmin + delta_R;
    
    dz = b.eta*dr;
    zmin = -(b.k2-1)*dz/2;  % assume zonull = 0
    zmax = -zmin;
    
    disp(['Grid RZ domain [Rmin,Rmax,Zmin,Zmax] = ',num2str([rmin,rmax,zmin,zmax])])
    
    % info for full field period
    phimin_fp = 0;
    phimax_fp = 2*pi/b.nperio;
    delta_phi_fp = phimax_fp - phimin_fp;
    dphi_fp = delta_phi_fp/b.ialfa;
    
    % for half field period
    phimin_hfp = 0;
    phimax_hfp = pi/b.nperio;
    delta_phi_hfp = phimax_hfp - phimin_hfp;
    dphi_hfp = delta_phi_hfp/(b.iald21 - 1);
    
    % build hfp grid
    nphi_hfp = b.iald21;
    nr = b.k2;
    nz = b.k2;
    
    phi_array = linspace(phimin_hfp,phimax_hfp,nphi_hfp);
    r_array = linspace(rmin,rmax,nr);
    z_array = linspace(zmin,zmax,nz);
end

idiv = 0;

ng_rz = b.k2;       % number of grid points per >>>>>>>>>> half FP <<<<<<<<
ng_phi = b.iald21;
nfp = b.nperio;
ig_ax = b.knull + 1;

% number of cells
nc_phi = b.ialfa;  % for full period

dphi_period = 2*pi/nfp;

% Set grid step sizes
dphi = dphi_period/nc_phi;
dr     = b.ronull/b.knull;
dz     = dr*b.eta;
dx    = 0.2*dr;  % step size for points around evaluation point +/-

%  Convert evaluation point to cartesian coordinates
cosphi0 = cos(xvec(2));
sinphi0 = sin(xvec(2));
xeval = xvec(1)*cosphi0;
yeval = xvec(1)*sinphi0;
zeval = xvec(3);

% 
%  Set grid points
% 
xpt(1) = xeval + dx; ypt(1) = yeval;      zpt(1) = zeval;
xpt(2) = xeval - dx; ypt(2) = yeval;      zpt(2) = zeval;
xpt(3) = xeval;      ypt(3) = yeval + dx; zpt(3) = zeval;
xpt(4) = xeval;      ypt(4) = yeval - dx; zpt(4) = zeval;
xpt(5) = xeval;      ypt(5) = yeval;      zpt(5) = zeval + dx;
xpt(6) = xeval;      ypt(6) = yeval;      zpt(6) = zeval - dx;

bxt = zeros(1,6);
byt = zeros(1,6);
bzt = zeros(1,6);

%
% loop over grid points and evaluate B
%
for i=1:6
  x = xpt(i);
  y = ypt(i);
  z = zpt(i);
  
  phi = atan2(y,x);
%   phi = phi + pi - pi*sign(phi); % phi must be in range 0:2pi now
%   nf = floor(phi/dphi_period); 
%   phi = phi - nf*dphi_period;  % map to field period phi
  phi = mod(phi,dphi_period);

  ray    = sqrt(x*x + y*y);
  cosphi  = x/ray;
  sinphi  = y/ray;

  af     = phi/dphi + 1; 
  lf     = floor(af);
  fq     = af - lf;
  lf1    = lf + 1;
  if lf1 > b.ialfa
      lf1 = 1;
  end
  ar     = (ray-b.rnull)/dr + ig_ax; 
  lr     = floor(ar);
  rq     = ar - lr;
  lr1    = lr + 1;
  as     = z/dz + ig_ax;
  ls     = floor(as);
  sq     = as - ls;
% *************************************************
  if (idiv == 1 || lr*(lr-ng_rz) >= 0 || ls*(ls-ng_rz) >= 0)
    idiv = 1;
%         write(*,*) 'idiv set 1 --> outside region?'
  else
    v1     = 1;
    v2     = 1;
    ls1    = ls;
    ls11   = ls1 + 1;
    ls2    = ls1;
    ls21   = ls11;
    if (lf > ng_phi)
      lf     = b.ialfa - lf + 2;
      v1     = -1;
      ls1    = -ls1  + ng_rz + 1;
      ls11   = -ls11 + ng_rz + 1;
    end
    if (lf1 > ng_phi)
      lf1    = b.ialfa - lf1 + 2;
      v2     = -1;
      ls2    = -ls2  + ng_rz + 1;
      ls21   = -ls21 + ng_rz + 1;
    end
%     *************************************************
%       interpolation : bilinear interpolation in r and z
%                       and linear in phi for each point
%     *************************************************
    fq1   = 1 - fq;
    fqv2  = fq*v2;
    fq1v1 = fq1*v1;
    a1    = (1-rq)*(1-sq);
    a2    = (1-rq)*   sq;
    a3    =    rq *(1-sq);
    a4    =    rq *   sq;
%     *************************************************       
    bri =  a1*( fq1v1*b.brg(lf ,lr ,ls1 ) +fqv2*b.brg(lf1,lr ,ls2 ) ) ...
          +a2*( fq1v1*b.brg(lf ,lr ,ls11) +fqv2*b.brg(lf1,lr ,ls21) ) ...
          +a3*( fq1v1*b.brg(lf ,lr1,ls1 ) +fqv2*b.brg(lf1,lr1,ls2 ) ) ...
          +a4*( fq1v1*b.brg(lf ,lr1,ls11) +fqv2*b.brg(lf1,lr1,ls21) );

    bsi =  a1*( fq1  *b.bzg(lf ,lr ,ls1 ) +fq  *b.bzg(lf1,lr ,ls2 ) ) ...
          +a2*( fq1  *b.bzg(lf ,lr ,ls11) +fq  *b.bzg(lf1,lr ,ls21) ) ...
          +a3*( fq1  *b.bzg(lf ,lr1,ls1 ) +fq  *b.bzg(lf1,lr1,ls2 ) ) ...
          +a4*( fq1  *b.bzg(lf ,lr1,ls11) +fq  *b.bzg(lf1,lr1,ls21) );
    bfi =  a1*( fq1  *b.bfg(lf ,lr ,ls1 ) +fq  *b.bfg(lf1,lr ,ls2 ) ) ...
          +a2*( fq1  *b.bfg(lf ,lr ,ls11) +fq  *b.bfg(lf1,lr ,ls21) ) ...
          +a3*( fq1  *b.bfg(lf ,lr1,ls1 ) +fq  *b.bfg(lf1,lr1,ls2 ) ) ...
          +a4*( fq1  *b.bfg(lf ,lr1,ls11) +fq  *b.bfg(lf1,lr1,ls21) );
%     *************************************************
    bxt(i) = -bfi*sinphi +bri*cosphi;
    byt(i) = +bfi*cosphi +bri*sinphi;
    bzt(i) =  bsi;
  end
end%   i=1,6
% *************************************************

% 
%  Calculate interpolated bvector
% 
bx = mean(bxt);
by = mean(byt);
bz = mean(bzt);

%  Calculate bval

% *************************************************
if ( bx == 0 && by == 0 && bz == 0 )
  idiv = 1;
  bval(1:3) = NaN;
%   write(*,*) 'idiv set 2 --> bmod == 0'
else
% *************************************************
  bval(1)  =  bx*cosphi0 + by*sinphi0;    % br
  bval(2)  = -bx*sinphi0 + by*cosphi0;    % bphi
  bval(3)  =  bz;
end


