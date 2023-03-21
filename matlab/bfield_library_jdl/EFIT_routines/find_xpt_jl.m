function xpt_info = find_xpt_jl(g,second,refine,tol,quiet,rguess,zguess,de)
%xpt_info = find_xpt_jl(g,second,refine,tol,quiet,rguess,zguess,de)
% tol is magnitude of bp at xpt
if nargin < 5 
    quiet = 0;
end
if nargin < 4
    tol = 1e-6;
end
if nargin < 3
    refine = 0;
end
if nargin < 2
    second = 0;
end
if nargin < 7
    rguess = [];
    zguess = [];
end
if nargin < 8
    de = [];
end

if isempty(de)
    der = (g.r(end)-g.r(1))*.05;
    dez = (g.z(end)-g.z(1))*.05; 
else
    der = de;
    dez = de;
end

first_guess_Bp_tol = 1e-3;  % If guess from boundary is > this then try another method (for limited configs)
Rmin_eval = g.r(3) + 1e-3;
Zmin_eval = g.z(3) + 1e-3;
Rmax_eval = g.r(end-2) - 1e-3;
Zmax_eval = g.z(end-2) - 1e-3;

niter_max = 15;

der_save = der;
dez_save = dez;

if isempty(rguess)

    % ----> Make guess from boundary. Note this won't work if the configuration is limited
    
    % Throw away zero points
    if ~isempty(g.bdry)
        ik = find(g.bdry(1,:) > 1.d-4);
        rb = g.bdry(1,ik);
        zb = g.bdry(2,ik);
    else
        rb = g.rzero; % This is not a good guess
        zb = 0; 
    end
    
    % Get minimum poloidal field on boundary (initial guess)
    b=bfield_geq_bicub(g,rb,zb);
    bp = sqrt(b.br.^2 + b.bz.^2);
    [bpx,ix] = min(bp);
    
    if bpx > first_guess_Bp_tol 
        if quiet == 0
            fprintf('Guess from boundary exceeded tolerance!\n')
        end
        der = (g.r(end)-g.r(1))*.15;
        dez = (g.z(end)-g.z(1))*.1;
    end
    rx = rb(ix);
    zx = zb(ix);
    if quiet == 0
        disp(['  Very rough first X-point. [Bp,R,Z] = ',num2str([bpx,rx,zx])])
    end
else
 
    % ----> Initial guess is supplied
    rx = rguess;
    zx = zguess;
    b=bfield_geq_bicub(g,rx,zx);
    bpx = sqrt(b.br.^2 + b.bz.^2);    
    if quiet == 0
        disp(['  Very rough first X-point is supplied: [Bp,R,Z] = ',num2str([bpx,rx,zx])])
    end
end
    

if refine == 1
  err = bpx;
%   de1 = 0.3;
  n1 = 100;  % grid dimension (square)
  bp = zeros(n1,n1);
  rg = zeros(n1,n1);
  zg = zeros(n1,n1);
  niter = 1;
  while err > tol && niter < niter_max
    rt = linspace(max([Rmin_eval,rx-der]),min([Rmax_eval,rx+der]),n1);
    zt = linspace(max([Zmin_eval,zx-dez]),min([Zmax_eval,zx+dez]),n1);
    for i = 1:n1 
      ztmp = repmat(zt(i),1,n1);
      b = bfield_geq_bicub(g,rt,ztmp);
      rg(i,:) = rt;
      zg(i,:) = ztmp;
      bp(i,:) = sqrt(b.br.^2 + b.bz.^2);      
    end
    [~,ix] = min(bp);
    [bpx,jx] = min(min(bp));
    ix=ix(jx);
    
    der = abs(rg(ix,jx) - rx);
    dez = abs(zg(ix,jx) - zx);
    
    err = bpx;
           
    rx = rg(ix,jx);
    zx = zg(ix,jx);

    niter = niter + 1;    
    
    if niter >= niter_max
        warning('niter_max exceeded for 1st x-point')
    end

  end
  if quiet == 0 
    disp(['  First Xpoint. [Bp,R,Z] = ',num2str([bpx,rx,zx])]);
  end
end


% Find second xpoint
rx2=[];
zx2=[];
der = der_save;
dez = dez_save;

if second == 1 
    
  % Guess that config is ~up-down symmetric
  rx2 = rx;
  zx2 = -zx;
  b = bfield_geq_bicub(g,rx2,zx2);
  bpx2 = sqrt(b.br.^2+b.bz.^2);
  err = bpx2;
  
  if refine == 1      
%       de2 = 0.25;
      n1 = 100;  % grid dimension (square)
      bp = zeros(n1,n1);
      rg = zeros(n1,n1);
      zg = zeros(n1,n1);
      niter = 1;
      while err > tol && niter < niter_max
    rt = linspace(max([Rmin_eval,rx2-der]),min([Rmax_eval,rx2+der,n1]));
    zt = linspace(max([Zmin_eval,zx2-dez]),min([Zmax_eval,zx2+dez,n1]));
          for i = 1:n1
              ztmp = repmat(zt(i),1,n1);
              b = bfield_geq_bicub(g,rt,ztmp);
              rg(i,:) = rt;
              zg(i,:) = ztmp;
              bp(i,:) = sqrt(b.br.^2 + b.bz.^2);
          end
          [~,ix] = min(bp);
          [bpx2,jx] = min(min(bp));
          ix=ix(jx);

          der = abs(rg(ix,jx) - rx2);
          dez = abs(zg(ix,jx) - zx2);

          err = bpx2;
          
          rx2 = rg(ix,jx);
          zx2 = zg(ix,jx);
          niter = niter + 1;
          if niter >= niter_max
              warning('niter_max exceeded for 2nd x-point')
          end
      end
      if quiet == 0
          disp(['  Second Xpoint. [Bp,R,Z] = ',num2str([bpx2,rx2,zx2])]);
      end
  end


end

xpt_info.rx = rx;
xpt_info.zx = zx;
xpt_info.bpx = bpx;
if second
    xpt_info.rx2 = rx2;
    xpt_info.zx2 = zx2;
    xpt_info.bpx2 = bpx2;
end

