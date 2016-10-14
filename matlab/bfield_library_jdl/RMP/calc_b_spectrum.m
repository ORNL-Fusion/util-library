% function calc_b_spectrum
clearvars;
tic;
XPLOT = 1; % 1 = psiN, 2 = sqrt(psiN)
INTERP = 0; % if 1 interpolate Br_mn plot 
NARDON = 0;

mmax = 20;
pnwant = linspace(0.1,0.998,10); %[0.99,0.995];
ntheta = 10*mmax;
nn = 3;
nphi = 8*nn;

nsurf = length(pnwant);

%
%  SET UP BFIELD
% 
% taper = 2900*[-1  1 -1  1 -1  1 -1  1 -1  1 -1  1]; % Even
taper = 2900*[-1  1 -1  1 -1  1  1 -1  1 -1  1 -1]; % Odd
rmp = build_d3d_icoils_jl(taper);
% gfile_name = 'C:\Work\DIII-D\164723\g164723.03059_410';
gfile_name = 'C:\Work\DIII-D\160884\efits\g160884.05009_537';
g = readg_g3d(gfile_name);
bfield.type = 'just_coils';
bfield.g = g;
bfield.coil = rmp.coil;
bfield.current = rmp.current;



phi = 2*pi*(0:nphi-1).'/nphi;
theta = 2*pi*(0:ntheta-1).'/ntheta;
dphi = phi(2) - phi(1);
dtheta = theta(2) - theta(1);


marr = (0:2*mmax) - mmax;

% bnorm = zeros(ntheta,nphi);
br_c = zeros(nsurf,2*mmax+1);
br_s = zeros(nsurf,2*mmax+1);
S = zeros(nsurf,1);

dpsidr = zeros(nsurf,ntheta);
dpsidz = zeros(nsurf,ntheta);


% Calculate pest coords for each surface, then Br_mn
pest = get_pest_coords(g,pnwant,ntheta);
for i = 1:nsurf
    [~,dpsidr(i,:),dpsidz(i,:)] = get_psi_bicub(g,pest.r(i,:),pest.z(i,:));
%     psiN(i,:) = pnwant(i)*ones(1,pest.ntheta);
end
if NARDON == 1
%     psiN = zeros(nsurf,ntheta);
    psiN = repmat(pnwant,ntheta,1).';
    psin12 = sqrt(psiN);
    ds12dr = ((-dpsidr-g.ssimag)/(g.ssibry-g.ssimag))./(2*psin12);
    ds12dz = ((-dpsidz-g.ssimag)/(g.ssibry-g.ssimag))./(2*psin12);
    ds12ave = sqrt(ds12dr.^2 + ds12dz.^2);
end

rn = dpsidr./sqrt(dpsidr.^2 + dpsidz.^2); % Unit vector in grad(psi)
zn = dpsidz./sqrt(dpsidr.^2 + dpsidz.^2);


for i = 1:nsurf
    S(i) = sum(pest.jac(i,:))*2*pi*dtheta;
    for k = 1:nphi
        phiarr(1,1:ntheta) = phi(k);
        [Bout,ierr] = bfield_general_rzphi(pest.r(i,:),pest.z(i,:),phiarr,bfield);
        if NARDON == 1
            Bax = bfield_geq_bicub(g,pest.r(i,:),pest.z(i,:));
            bnorm{i}(:,k) = (Bout.br.*ds12dr(i,:).' + Bout.bz.*ds12dz(i,:).')./(Bax.bphi./pest.r(i,:).');
        else
            bnorm{i}(:,k) = Bout.br.*rn(i,:).' + Bout.bz.*zn(i,:).';
        end
    end
end

for i = 1:nsurf
    for mind = 1:2*mmax + 1
        for k = 1:nphi
            alpha_mn = nn*phi(k) - marr(mind)*theta;
            br_c(i,mind) = br_c(i,mind) + 2*dtheta*dphi*sum(pest.jac(i,:).'.*bnorm{i}(:,k).*cos(alpha_mn))./S(i); %A.15
            br_s(i,mind) = br_s(i,mind) + 2*dtheta*dphi*sum(pest.jac(i,:).'.*bnorm{i}(:,k).*sin(alpha_mn))./S(i);
        end
    end
end
if NARDON
    dsave = mean(ds12ave(i,:));
    br_c(i,:) = 2*(br_c(i,:))/(g.rzero*dsave);
    br_s(i,:) = 2*(br_s(i,:))/(g.rzero*dsave);
end
Br_mn = sqrt(br_c.^2 + br_s.^2); % After A.15

% Bres and island widths
mlow = ceil(nn*min(pest.q));
mhigh = floor(nn*max(pest.q));
mres= mlow:mhigh;

numres = length(mres);
if numres < 2 
    error('Extend this')
end
widths = zeros(numres,1);
pnres  = zeros(numres,1);
bres   = zeros(numres,1);
wid    = zeros(numres,1);

shear = interp1(g.pn,deriv(g.pn,g.qpsi),pest.psiN);
dpsi = abs(g.ssibry - g.ssimag);
for i = 1:numres
    qres = mres(i)/nn;
    pnres(i) = interp1(pest.q,pest.psiN,qres);
    sres = interp1(pest.psiN,shear,pnres(i));
    SAres = interp1(pest.psiN,S,pnres(i));
    mind = find(marr == mres(i));
    if length(mind) ~= 1 
        error('m not found')
    end
    bres(i) = interp1(pest.psiN,Br_mn(:,mind),pnres(i));
    wid(i) = sqrt((16/(mres(i)*dpsi))*(qres/sres)*(SAres/(4*pi*pi))*bres(i));
end
pnchir = 0.5*(pnres(2:numres) + pnres(1:numres-1));
dpchir = pnres(2:end) - pnres(1:end-1);
dwchir = 0.5*(wid(2:end) + wid(1:end-1));
chir = dwchir./dpchir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   MAKE PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if XPLOT == 1
    radarr = pnwant;
    radlab = '\psi_N';
elseif XPLOT == 2
    radarr = sqrt(pnwant);
    radlab = '\psi_N^{1/2}';
end

fontsize = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  ISLAND WIDTHS ON q PROFILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on; box on;
plot(pest.psiN,pest.q,'b-','linewidth',2)
for i = 1:numres
    qtmp = interp1(pest.psiN,pest.q,pnres(i));
    plot([pnres(i)-0.5*wid(i),pnres(i)+0.5*wid(i)],[qtmp,qtmp],'k-','linewidth',3)
end
xlabel(radlab,'fontsize',fontsize)
ylabel('q','fontsize',fontsize)
set(gca,'fontsize',fontsize)
title(sprintf('n=%d island widths',nn))
axis([0,1,0,ceil(max(pest.q))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Chirikov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on; box on;
plot(pnchir,chir,'ko')
plot([0,1],[1,1],'k-')
xlabel(radlab,'fontsize',fontsize)
ylabel('\sigma_{chir}','fontsize',fontsize)
set(gca,'fontsize',fontsize)
title(sprintf('n=%d Chirikov Parameter',nn))
axis([0,1,0,ceil(max(chir))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Br_mn contour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on; box on;
if INTERP
    h = pcolor(marr,radarr,Br_mn*1e4);
else
    h = pcolor(marr-.5,radarr,Br_mn*1e4);
end
% h = surf(marr,xarr,br*1e4);

set(h,'edgecolor','none');
if INTERP
    set(h,'facecolor','interp');
end
plot(-pest.q*nn,radarr,'w--','linewidth',3)
xlabel('Poloidal mode number m','fontsize',fontsize)

ylabel(radlab,'fontsize',fontsize)
colorbar;
title(sprintf('B_r^{n=%d} [G]',nn))
% axis([-20,0,0.4,1])
% axis([0,20,0.4,1])
colormap(colorflipper(1024,'jet_to_white'))
% colormap(colorflipper(1024,'jet'))

% % OVERSAMPLE
% mtest = linspace(-20,0.1,2000);
% stest = linspace(0.7,0.998,4000);
% s2test = stest.^2;
% [mm,ss] = meshgrid(mtest,stest);
% 
% brsamp = interp2(marr,sqrt(pnwant),br,mm,ss);
% 
% figure; hold on; box on;
% h = pcolor(mtest,stest,brsamp*1e4);
% set(h,'edgecolor','none','facecolor','interp');
% plot(-g.qpsi*nn,sqrt(g.pn),'w')
% xlabel('Poloidal mode number m','fontsize',14)
% set(gca,'fontsize',14)
% ylabel('\psi_N^{1/2}','fontsize',14)
% colorbar;
% title(sprintf('B_r_{(m,%d)} [Gauss]',nn))
% axis([-20,0,0.7,1])
toc
