function [s,u,v,dRds,dZds,dRdu,dRdv,dZdu,dZdv]=rzp_to_suv(R,Z,phi,wout,tolfun,tolx,quiet)
if nargin < 5 | isempty(tolfun)
    tolfun = 1e-10;
end
if nargin < 6 | isempty(tolx)
    tolx = 1e-10;
end
if nargin < 7
   quiet = 1;
end
sstep = 1e-8;

% R = 4;
% Z = 0.1;
% phi = 0.1;
% wout_file = 'C:\Work\Pellets\wout_vmec_pcnet.nc';
% wout=load_wout(wout_file);
npts = length(R);
if quiet 
    qval = 'off';
else
    qval = 'final';
end
opts=optimoptions('lsqnonlin','TolFun',tolfun,'TolX',tolx,'Display',qval);
opts.TolFun=tolfun;
opts.TolX=tolx;

% x0 = [0.5,0.1];
x0=repmat([0.5,0.1],npts,1);
xub=repmat([1,2*pi],npts,1);   %s,u
xlb=repmat([0,-2*pi],npts,1);  %s,u
% % xub = [1,2*pi];
% % xlb = [0,-2*pi];
xfinal=lsqnonlin(@rzp_to_suv_minfun,x0,xlb,xub,opts);
s = xfinal(:,1).';
u = xfinal(:,2).';
v=phi;

if nargout > 3
    [R,Z,dRdu,dRdv,dZdu,dZdv] = suv_to_rz(s,u,phi,wout);
    [R2,Z2] = suv_to_rz(s+sstep,u,phi,wout);
    dRds = (R2-R)/sstep;
    dZds = (Z2-Z)/sstep;
end

    function f=rzp_to_suv_minfun(x)        
        stmp=x(:,1);
        utmp=x(:,2);
        [Rout,Zout] = suv_to_rz(stmp.',utmp.',phi,wout);        
        f = [Rout-R,Zout-Z];        
    end
end