% function calc_b_spectrum
clearvars;
tic;
XPLOT = 1; % 1 = psiN, 2 = sqrt(psiN)
INTERP = 0;

mmax = 32;
pnwant = linspace(0.1,0.998,256); %[0.99,0.995];
ntheta = 100*mmax;
nn = 3;
nphi = 20*nn;

nsurf = length(pnwant);

%
%  SET UP BFIELD
% 
taper = 2900*[-1  1 -1  1 -1  1 -1  1 -1  1 -1  1]; % Even
% taper = 2900*[-1  1 -1  1 -1  1  1 -1  1 -1  1 -1]; % Odd
rmp = build_d3d_icoils_jl(taper);
gfile_name = 'C:\Work\DIII-D\164723\g164723.03059_410';
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
psiN = zeros(nsurf,ntheta);

% Calculate pest coords for each surface
pest = get_pest_coords(g,pnwant,ntheta);
for i = 1:nsurf
    [~,dpsidr(i,:),dpsidz(i,:)] = get_psi_bicub(g,pest.r(i,:),pest.z(i,:));
    psiN(i,:) = pnwant(i)*ones(1,pest.ntheta);
end
rn = dpsidr./sqrt(dpsidr.^2 + dpsidz.^2); % Unit vector in grad(psi)
zn = dpsidz./sqrt(dpsidr.^2 + dpsidz.^2);

for i = 1:nsurf
    S(i) = sum(pest.jac(i,:))*2*pi*dtheta;
    for k = 1:nphi
        phiarr(1,1:ntheta) = phi(k);
        [Bout,ierr] = bfield_general_rzphi(pest.r(i,:),pest.z(i,:),phiarr,bfield);
        bnorm{i}(:,k) = Bout.br.*rn(i,:).' + Bout.bz.*zn(i,:).';
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

Br_mn = sqrt(br_c.^2 + br_s.^2); % After A.15

if XPLOT == 1
    yarr = pnwant;
    ylab = '\psi_N';
elseif XPLOT == 2
    yarr = sqrt(pnwant);
    ylab = '\psi_N^{1/2}';
end

figure; hold on; box on;
if INTERP
h = pcolor(marr,yarr,Br_mn*1e4);
else
    h = pcolor(marr-.5,yarr,Br_mn*1e4);
end
% h = surf(marr,xarr,br*1e4);

set(h,'edgecolor','none');
if INTERP
    set(h,'facecolor','interp');
end
plot(-pest.q*nn,yarr,'w--','linewidth',3)
xlabel('Poloidal mode number m','fontsize',14)
set(gca,'fontsize',14)
ylabel(ylab,'fontsize',14)
colorbar;
title(sprintf('B_r^{n=%d} [G]',nn))
axis([-20,0,0.4,1])
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
