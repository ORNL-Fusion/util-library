clearvars;

%% Options
PLOT_GFILE = 0;
RADIAL_METHOD = 'tanh';  % Options: 'psiN', 'tanh'
POLOIDAL_METHOD = 'theta'; % Options: 'uniform', 'theta', 'parallel' (must use maskLim == 1)
maskPFR = 1;           % If 1 then treat PFR psiN as 2 - psiN
maskPoloidalBySep = 1; % If 1 then only apply poloidal function for psiN >= 1
maskLim = 1;           % If 1 then set fun = 0 outside of g.lim

%% Settings
gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\PSIscidac\WEST\g055866.47179_V2';  % efit file

nR = 30*2;  % Grid resolution to build surface
nZ = 40*2;

cTanh_ne = [1 0.04 0.3 0.1 3e-3 3e-5 0 -0.01];  % tanh coefficents, only used if RADIAL_METHOD == 'psiN'

Apol = 0.95; % theta method is applied a 1 + Apol*sin(theta)

Ztest = linspace(-0.6,0.6,100); % Line over which to evaluate 2D fit
Rtest = 2.95*ones(size(Ztest));

%% Load and optionally plot gfile, find xpt (for PFR filtering)
g = readg_g3d(gfile_name);
Xpt = find_xpt_jl(g,1,1,1e-10,1);
if PLOT_GFILE
    plot_gfile(g);
end

%% Set up interpolation mesh
R1D = linspace(g.r(3),g.r(end-2),nR);
Z1D = linspace(g.z(3),g.z(end-2),nZ);
[R2D,Z2D] = meshgrid(R1D,Z1D);

%% Compute some quantities
psiN = reshape(calc_psiN(g,R2D(:),Z2D(:)),nZ,nR);
theta = atan2(Z2D - g.zmaxis,R2D - g.rmaxis);

%% Adjust PFR
if maskPFR
    these = psiN < 1 & Z2D < Xpt.zx;
    psiN(these) = 2 - psiN(these);
end

%% Choose radial fun
switch RADIAL_METHOD
    case 'psiN'
        radialFun_ne = psiN;
        % put your psiN to temperature/density mapping here
    case 'tanh'
        radialFun_ne = evaluate_tanh_fit(cTanh_ne,psiN);
    otherwise
        error('Unknown RADIAL_METHOD: %s\n',RADIAL_METHOD)
end

%% Choose poloidal fun
switch POLOIDAL_METHOD
    case 'uniform'
        poloidalFun_ne = ones(nZ,nR);
    case 'theta'
        poloidalFun_ne = (1 + Apol*sin(theta));
        poloidalFun_ne = poloidalFun_ne./(max(poloidalFun_ne(:)));
    case 'parallel'
        bfield.type = 'gfile';
        bfield.nsym = 1;
        bfield.g = g;
        Lmax = 100;
        dL = 0.01;
        for i = 1:nR
            for j = 1:nZ
                s = follow_fieldlines_rzphi_dl(bfield,R2D(j,i),Z2D(j,i),0,dL,Lmax/dL);
                dsfsdaf
                aaa=1
            end
        end
    otherwise
        error('Unknown POLOIDAL_METHOD: %s\n',POLOIDAL_METHOD)
end

%% Apply only outside separatrix
if maskPoloidalBySep
    these = psiN < 1;
    poloidalFun_ne(these) = 1;
end

%% Define total function
myFun = radialFun_ne.*poloidalFun_ne;

if maskLim
    mask = ~inpolygon(R2D,Z2D,g.lim(1,:),g.lim(2,:));
    myFun(mask) = 0;
end

%% Make cells
for j = 1:nR - 1
    for i = 1:nZ - 1
        rCell(:,(i-1)*(nZ-1) + j) = [R2D(i,j),R2D(i+1,j),R2D(i+1,j+1),R2D(i,j+1)];
        zCell(:,(i-1)*(nZ-1) + j) = [Z2D(i,j),Z2D(i+1,j),Z2D(i+1,j+1),Z2D(i,j+1)];
        vCell((i-1)*(nZ-1) + j,1) = myFun(i,j);
    end
end

%% Plot contours
figure;
hold on; box on; set(gcf,'color','w');set(gca,'fontsize',14); grid on;
axis equal;
patch(rCell,zCell,vCell,'edgecolor','none')
% h = pcolor(R2D(~mask),Z2D(mask),myFun(mask));
% set(h,'edgecolor','none');

colorbar;
c = parula;
c(1,:) = [1 1 1];
colormap(c);
plot(g.lim(1,:),g.lim(2,:),'k')
plot(g.bdry(1,:),g.bdry(2,:),'k')


%% interpolate onto line
if 0

    plot(Rtest,Ztest,'r','linew',3)
    funTest = interp2(R2D,Z2D,myFun,Rtest,Ztest,'makima');

    figure; hold on; box on; set(gcf,'color','w');set(gca,'fontsize',14); grid on;
    plot(Ztest,funTest)
    xlabel('Z (m)')
    ylabel('test function')
end