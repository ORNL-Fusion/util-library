clearvars;

% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\PSIscidac\WEST\g055866.47179_V2';  % efit file
gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\NegD\current scan\194288_3000_1.gfile';
g = readg_g3d(gfile_name);

%% Make rectangular grid to refine triangles
nR = 30*3;  % Grid resolution to build surface
nZ = 40*3;

R1D = linspace(g.r(3),g.r(end-2),nR);
Z1D = linspace(g.z(3),g.z(end-2),nZ);
[R2D,Z2D] = meshgrid(R1D,Z1D);

% Toss points outside of lim
mask = inpolygon(R2D,Z2D,g.lim(1,:),g.lim(2,:));
R2D = R2D(mask);
Z2D = Z2D(mask);

%% Refine limiter
P = g.lim';
P = refineLim(P,3e-2);

% and get rid of repeated points
tolMin = 1e-6;
dx = sqrt(diff(P(:,1)).^2 + diff(P(:,2)).^2 );
dx(end+1) = sqrt( (P(1,1) - P(end,1))^2 + (P(1,2) - P(end,2))^2 );
P(dx < tolMin,:) = [];

% Constrain edges
n = size(P,1);
for i = 1:n
    if i == n
        C(i,:) = [i,1];
    else
        C(i,:) = [i,i+1];
    end
end
% C = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 1];

%% Toss it all in the blender
P = [P;[R2D,Z2D]];

DT = delaunayTriangulation(P,C);
IO = isInterior(DT);
TR = triangulation(DT.ConnectivityList(IO,:),DT.Points);

figure; hold on; box on; grid on; set(gcf,'color','w');

triplot(TR)

plot(g.lim(1,:),g.lim(2,:),'k','linew',3)
% plot(R2D(:),Z2D(:),'.')


%% Make cells
conn = DT.ConnectivityList(IO,:);
nTri = size(conn,1);
for i = 1:nTri    
    RTri(:,i) = [DT.Points(conn(i,1),1),DT.Points(conn(i,2),1),DT.Points(conn(i,3),1)];
    ZTri(:,i) = [DT.Points(conn(i,1),2),DT.Points(conn(i,2),2),DT.Points(conn(i,3),2)];
end


%% Get centers
Xpt = find_xpt_jl(g,1,1,1e-10);
Centers = incenter(TR);
psiN = calc_psiN(g,Centers(:,1),Centers(:,2));
maskPFR = 1;
if maskPFR
    these = psiN < 1 & Centers(:,2)' < Xpt.zx;
    psiN(these) = 2 - psiN(these);
end

% cTanh_ne = [1 0.04 0.3 0.1 3e-3 3e-5 0 -0.01];  % tanh coefficents, only used if RADIAL_METHOD == 'psiN'
% cTanh_te = [1 0.04 0.3 0.1 3e-3 3e-5 0 -0.01];  % tanh coefficents, only used if RADIAL_METHOD == 'psiN'
cTanh_ne = [0.9842    0.0336    0.6436    0.1551   -0.0044    0.0004   -0.0000   -0.001   -0.00];
% cTanh_te = [0.9641    0.0458    0.5112    0.0198    0.3138   -0.0076    0.0002   -0.0010    0.0008];
cTanh_te = [0.9606    0.0521    0.5740    0.0157    0.2345];

radialFun_ne = evaluate_tanh_fit(cTanh_ne,psiN)*1e20;
radialFun_te = evaluate_tanh_fit(cTanh_te,psiN)*1e3;

te = radialFun_te;
ne = radialFun_ne;

max(psiN)

pc = phys_const;
fname = 'C:\Users\jjl\Dropbox (ORNL)\ADAS\adf11_all\scd12\scd12_h.dat';
indRC = 1;
RateCoeff = read_adas_adf11_file(fname); 
sigma_v = 1e-6*interp_adas_rate_coefficient(te,ne./1e6,RateCoeff.te,RateCoeff.ne,squeeze(RateCoeff.coeff(:,:,indRC)) );
am = 2;
vp = sqrt(2.*3*pc.eV/(am*pc.mp));
%  
% plotThis(plotThis > 1e4) = 1e4; 
% plotThis = log10(plotThis); 
% PlotThis.myTitle = 'log10(\lambda_{mfp}) (m)';


figure; hold on; box on; grid on; set(gcf,'color','w');
plotThis = te;
patch(RTri,ZTri,plotThis,'edgecolor','none');
plot(g.lim(1,:),g.lim(2,:),'k','linew',3)
colorbar;
title('T_e eV')

figure; hold on; box on; grid on; set(gcf,'color','w');
plotThis = log10(ne);
patch(RTri,ZTri,plotThis,'edgecolor','none');
plot(g.lim(1,:),g.lim(2,:),'k','linew',3)
colorbar;
title('log10(n_e)')

figure; hold on; box on; grid on; set(gcf,'color','w');
% plotThis = nan(size(RTri,2),1);
% plotThis = te;
% plotThis = vp./(ne.*sigma_v);
plotThis = log10(vp./(ne.*sigma_v));
patch(RTri,ZTri,plotThis,'edgecolor','none');
plot(g.lim(1,:),g.lim(2,:),'k','linew',3)
colorbar;
title('log10(\lambda_{mfp})')



function P = refineLim(P,tol)
%% refine to meet some maximum edge length
% tolMax = 10e-2;
dx = sqrt(diff(P(:,1)).^2 + diff(P(:,2)).^2 );
Pnew(1,:) = P(1,:);
for i = 2:size(P,1)
    dxThis = dx(i-1);
    if dxThis <= tol
        ind = size(Pnew,1) + 1;
        Pnew(ind,1) = P(i,1);
        Pnew(ind,2) = P(i,2);
    else
        Rnew = linspace(P(i-1,1),P(i,1),ceil(dxThis/tol));
        Znew = linspace(P(i-1,2),P(i,2),ceil(dxThis/tol));

        ind = size(Pnew,1)+1:size(Pnew,1)+length(Rnew)-1;
        Pnew(ind,1) = Rnew(2:end);
        Pnew(ind,2) = Znew(2:end);
    end
end

P = Pnew;
end
