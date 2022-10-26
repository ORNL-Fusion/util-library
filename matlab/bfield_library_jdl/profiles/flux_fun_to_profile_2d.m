clearvars;

%% Options
PLOT_GFILE = 0;
RADIAL_METHOD = 'tanh';  % Options: 'psiN', 'tanh'
POLOIDAL_METHOD = 'parallel'; % Options: 'uniform', 'theta', 'parallel' (must use maskLim == 1)
maskPFR = 1;           % If 1 then treat PFR psiN as 2 - psiN
maskPoloidalBySep = 1; % If 1 then only apply poloidal function for psiN >= 1
maskLim = 1;           % If 1 then set fun = 0 outside of g.lim

%% Settings
gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\PSIscidac\WEST\g055866.47179_V2';  % efit file

nR = 30;  % Grid resolution to build surface
nZ = 40;

cTanh = [1 0.04 0.3 0.1 3e-3 3e-5 0 -0.01];  % tanh coefficents, only used if RADIAL_METHOD == 'psiN'

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
        radialFun = psiN;
        % put your psiN to temperature/density mapping here
    case 'tanh'
        radialFun = evaluate_tanh_fit(cTanh,psiN);
    otherwise 
        error('Unknown RADIAL_METHOD: %s\n',RADIAL_METHOD)
end

%% Choose poloidal fun 
switch POLOIDAL_METHOD
    case 'uniform'
        poloidalFun = ones(nZ,nR);
    case 'theta'
        poloidalFun = (1 + Apol*sin(theta));
        poloidalFun = poloidalFun./(max(poloidalFun(:)));
    case 'parallel'
        bfield.type = 'gfile';
        bfield.nsym = 1;
        bfield.g = g;
        Lmax = 100;
        dL = 0.01;
        for i = 1:nR
            for j = 1:nZ
                s = follow_fieldlines_rzphi_dl(bfield,R2D(j,i),Z2D(j,i),0,dL,Lmax/dL);
        aaa=1
            end
        end
    otherwise 
        error('Unknown POLOIDAL_METHOD: %s\n',POLOIDAL_METHOD)
end

%% Apply only outside separatrix
if maskPoloidalBySep
    these = psiN < 1;
    poloidalFun(these) = 1;
end

%% Define total function
myFun = radialFun.*poloidalFun;

if maskLim
    mask = ~inpolygon(R2D,Z2D,g.lim(1,:),g.lim(2,:));
    myFun(mask) = 0;
end

%% Plot contours
figure; hold on; box on; set(gcf,'color','w');set(gca,'fontsize',14); grid on; 
axis equal;
h = pcolor(R2D,Z2D,myFun);
set(h,'edgecolor','none');
colorbar;
c = parula;
c(1,:) = [1 1 1];
colormap(c);
plot(g.lim(1,:),g.lim(2,:),'k')
plot(g.bdry(1,:),g.bdry(2,:),'k')
plot(Rtest,Ztest,'r','linew',3)

%% interpolate onto line
funTest = interp2(R2D,Z2D,myFun,Rtest,Ztest,'makima');

figure; hold on; box on; set(gcf,'color','w');set(gca,'fontsize',14); grid on;
plot(Ztest,funTest)
xlabel('Z (m)')
ylabel('test function')