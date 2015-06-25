function xpt_info = find_xpt_jl(g,second,refine,tol,quiet,rguess,zguess,de)
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
    %   de1 = 0.3;
    %   de2 = 0.25;
    de1 = [];
    de2 = [];
else
    de1 = de;
    de2 = de;
end

if isempty(de1)
    de1 = (g.r(end)-g.r(1))*.15;
    de2 = (g.z(end)-g.z(1))*.15; 
end

niter_max = 15;

if isempty(rguess)

    % ----> Make guess from boundary. Note this won't work if the configuration is limited
    
    % Throw away zero points
    ik = find(g.bdry(1,:) > 1.d-4);
    rb = g.bdry(1,ik);
    zb = g.bdry(2,ik);
    
    % Get minimum poloidal field on boundary (initial guess)
    b=bfield_geq_bicub(g,rb,zb);
    bp = sqrt(b.br.^2 + b.bz.^2);
    [bpx,ix] = min(bp);
    rx = rb(ix);
    zx = zb(ix);
    if quiet == 0
        disp([' Very rough first X-point. [Bp,R,Z] = ',num2str([bpx,rx,zx])])
    end
else
 
    % ----> Initial guess is supplied
    rx = rguess;
    zx = zguess;
    b=bfield_geq_bicub(g,rx,zx);
    bpx = sqrt(b.br.^2 + b.bz.^2);    
    if quiet == 0
        disp([' Very rough first X-point is supplied: [Bp,R,Z] = ',num2str([bpx,rx,zx])])
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
    rt = linspace(rx-0.5*de1,rx+0.5*de1,n1);
    zt = linspace(zx-0.5*de1,zx+0.5*de1,n1);
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
    de1 = sqrt( (rg(ix,jx)-rx)^2 + (zg(ix,jx)-zx)^2);
    err = bpx;
    
    rx = rg(ix,jx);
    zx = zg(ix,jx);
    niter = niter + 1;    
    if niter >= niter_max
        warning('niter_max exceeded for 1st x-point')
    end

  end
  if quiet == 0 
    disp([' First Xpoint. [Bp,R,Z] = ',num2str([bpx,rx,zx])]);
  end
end


% Find second xpoint
rx2=[];
zx2=[];

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
          rt = linspace(rx2-0.5*de2,rx2+0.5*de2,n1);
          zt = linspace(zx2-0.5*de2,zx2+0.5*de2,n1);
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
          de2 = sqrt( (rg(ix,jx)-rx2)^2 + (zg(ix,jx)-zx2)^2);          
          err = bpx2;
          
          rx2 = rg(ix,jx);
          zx2 = zg(ix,jx);
          niter = niter + 1;
          if niter >= niter_max
              warning('niter_max exceeded for 2nd x-point')
          end
      end
      if quiet == 0
          disp([' Second Xpoint. [Bp,R,Z] = ',num2str([bpx2,rx2,zx2])]);
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

