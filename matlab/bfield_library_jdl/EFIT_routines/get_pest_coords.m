function pest = get_pest_coords(g,pnwant,ntheta,close)
if nargin < 1
    error('Must provide g file')
end
if nargin < 2
    pnwant = g.pn;
    pnwant(end) = 0.999;
end
if nargin < 3
    ntheta = 128;
end
if nargin < 4
    close = 0;
end

% clearvars;
% close = 0; % if 1 then phi goes from theta = [0,2*pi], else [0,2*pi)
% pnwant = [0.5,0.6,0.99];
% psiedge = 0.999;
% gfile_name = 'C:\Work\DIII-D\164723\g164723.03059_410';
% g = readg_g3d(gfile_name);


if any(pnwant > 1. | pnwant < 0.)
    error('PsiN must be < 1.')
end
npsin = length(pnwant);


rpest = zeros(npsin,ntheta);
zpest = zeros(npsin,ntheta);
jac = zeros(npsin,ntheta);

% Set up bfield and follow field lines
bfield.type = 'gfile';
bfield.g = g;

qpsi = interp1(g.pn,g.qpsi,pnwant);

[rstart,zstart] = calc_RZ_at_psiN_theta(g,pnwant,0.);
phistart = 0.;
dphi = 0.5*pi/180;
nsteps = round(1.1*2*pi/dphi*max(qpsi));

f = follow_fieldlines_rzphi_dphi(bfield,rstart,zstart,phistart,dphi,nsteps);
f.theta = atan2(f.z - g.zmaxis,f.r - g.rmaxis);
f.theta(f.theta<0) = f.theta(f.theta<0)+2*pi;
% max(f.theta)

theta_pest = zeros(npsin,ntheta);
for i = 1:npsin
    % Get one period in theta: 0:2*pi
    if f.theta(2,i) < f.theta(1,i) 
        error('Need to take this case into account')
    end
    icross = find(f.theta(2:end,i) - f.theta(1:end-1,i) < 0.);
    if isempty(icross)
        error('Probably need to increase nsteps')
    end
    th_fl  = [f.theta(1:icross(1),i);f.theta(icross(1)+1,i)+2*pi];
%     figure; plot(th_fl); title((num2str(i)))
    phi_fl = f.phi(1:icross(1)+1);
%       figure; plot(phi_fl,th_fl); title((num2str(i)))
    phi_th0 = interp1(th_fl,phi_fl,2*pi); % interpolate to theta = 2*pi 
    r_fl = f.r(1:icross(1)+1,i);
    z_fl = f.z(1:icross(1)+1,i);   
    if close == 1
        theta_pest(i,:) = 2*pi*(0:ntheta-1)./(ntheta-1); 
    else
        theta_pest(i,:) = 2*pi*(0:ntheta-1)./ntheta;        
    end
    phi_pest = phi_th0*theta_pest(i,:)/2/pi;       
    rpest(i,:) = interp1(phi_fl,r_fl,phi_pest);
    zpest(i,:) = interp1(phi_fl,z_fl,phi_pest);
    
    b = bfield_geq_bicub(g,rpest(i,:),zpest(i,:));
    bpol = sqrt(b.br.^2 + b.bz.^2).';
    jac(i,:) = abs(qpsi(i)/(g.rzero*g.bcentr))*bpol.*rpest(i,:).^3;

end

pest.psiN = pnwant;
pest.r = rpest;
pest.z = zpest;
pest.jac = jac;
pest.q = qpsi;
pest.npsi = npsin;
pest.ntheta = ntheta;
pest.theta = theta_pest;
