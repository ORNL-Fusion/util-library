function poinc = make_poincare_coils_only(coil,Rstart,Zstart,phistart,phi_want,npoints_want,dphi,sort_it)
% R,Z = [m], phi = [radians]
% R,Z can be arrays
if nargin == 0
    phistart = 0;
    Rstart = 6;
    Zstart = 0;   
end

if nargin < 5
    phi_want = phistart;
end
if nargin < 6
    npoints_want = 10;
end
if nargin < 7
    dphi = 0.5*pi/180;
end
if nargin < 8
    sort_it = 0;
end


tic;
fprintf('Making %d poincare surfaces\n',length(Rstart))

%--------------------------------------------------------------------------
%      BFIELD SETUP
%--------------------------------------------------------------------------

if nargin == 0    
    coils_file = 'C:\Work\Stellarator\W7X EMC3 modeling\Mark coils and input\coils.w7x';            
    coil = load_vmec_coils_file(coils_file);    
    taper = [13067,12931,13203,14564,14700,8983,-3267,-2756,2756];  % 0kA mimic    
    coil = set_w7x_current(coil,taper); % taper = [I1,I2,I3,I4,I5,IA,IB,IS1,IS2];
    fprintf('Bfield setup took %f seconds.\n',toc); tic;
end

bfield.type = 'just_coils';
bfield.coil = coil.coil;
bfield.current = coil.current;
bfield.nsym = coil.num_periods;

bounds.type = 'box';
bounds.Rlims = [5,6.7];
bounds.Zlims = [-1.5,1.5];
bfield.bounds = bounds;

%--------------------------------------------------------------------------
%      POINCARE SETUP
%--------------------------------------------------------------------------


nsurfs = length(Rstart);
if nsurfs ~= length(Zstart)
    if length(Zstart) == 1
        Zstart = Zstart*ones(size(Rstart));
    else
        error('Input inconsistency')
    end
end
if length(phistart) ~= 1
    error('No')
end

NFP = bfield.nsym;
ntran = npoints_want/NFP;
nsteps_real = ntran*2*pi/dphi;
nsteps = round(nsteps_real);
if abs(nsteps_real - nsteps) > 1e-6
    error('Make this divisible')
end

nowarn = 0;
istart = (phi_want - phistart)/dphi + 1;
istride = 2*pi/NFP/dphi;
npoinc = length(istart:istride:nsteps+1);
Rpoinc = NaN(npoinc,nsurfs);
Zpoinc = NaN(npoinc,nsurfs);
for i = 1:nsurfs
    fprintf('Working on surface %d of %d\n',i,nsurfs);
    f = follow_fieldlines_rzphi_dphi(bfield,Rstart(i),Zstart(i),phistart,dphi,nsteps,nowarn);
    Rpoinc(:,i) = f.r(istart:istride:end,1);
    Zpoinc(:,i) = f.z(istart:istride:end,1);   
    if 0
        x = f.r.*cos(f.phi);
        y = f.r.*sin(f.phi);
        z = f.z;
        figure; hold on; box on;
        plot3(x,y,z)
    end
end


% Rpoinc = f.r(istart:istride:end,:);
% Zpoinc = f.z(istart:istride:end,:);
fprintf('Field line following took %f seconds.\n',toc); tic;

%--------------------------------------------------------------------------
%      PLOTTING
%--------------------------------------------------------------------------


if sort_it
%     phistart = 0;
    [Rax,Zax] = find_axis(bfield,phistart);
%     Rax = mean(mean(Rpoinc));
%     Zax = mean(mean(Zpoinc));
    tpoinc = atan2(Zpoinc-Zax,Rpoinc-Rax);
    rpoinc = sqrt((Rpoinc-Rax).^2 + (Zpoinc-Zax).^2);
    [tpoinc,Isort]=sort(tpoinc);                       %sort by increasing theta
    for j = 1:nsurfs
        rpoinc(:,j) = rpoinc(Isort(:,j),j);
        Rpoinc(:,j) = Rpoinc(Isort(:,j),j);
        Zpoinc(:,j) = Zpoinc(Isort(:,j),j);
    end
%     rpoinc = rpoinc(Isort);
%     Zpoinc = Zpoinc(Isort);
%     Rpoinc = Rpoinc(Isort);
end
    

figure; hold on; box on;
if sort_it    
    plot(Rpoinc,Zpoinc,'.-','linewidth',2)
else
    plot(Rpoinc,Zpoinc,'o')
end
xlabel('R [m]','fontsize',14)
ylabel('Z [m]','fontsize',14)
set(gca,'fontsize',14);
title(strcat('\Phi = ',sprintf('%f',phi_want*180/pi)));

poinc.Rpoinc = Rpoinc;
poinc.Zpoinc = Zpoinc;
poinc.phi = phi_want;