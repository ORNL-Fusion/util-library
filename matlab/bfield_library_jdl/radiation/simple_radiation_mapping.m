clearvars;
METHOD = 2;

tic;
PLOT_TYPE = 3; % 1 is power, 2 is power density (W/m), 3 power density (W/m^2)

%%
%--------------------------------------------------------------------------
% Define PFC and Source
% Note: PFC must be non-overlapping and simply sorted
%--------------------------------------------------------------------------

CASE = 3; % 1 = tokamak structure from gfile, 2 = simple defined geometry

switch CASE
    case 1
        gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Haskey\g179477.02280';
        g = readg_g3d(gfile_name);
        % remove weird point specific to this gfile
        g.lim(:,7) = [];
        pfc.R = g.lim(1,:);
        pfc.Z = g.lim(2,:);
% Simple routine to cleanup PFC (checks for duplicate points and closes
% contour)
pfc = cleanupPFC(pfc);
        Source.R = [1.5,1.45,2.3,2.3,2.3,1.5];
        Source.Z = [0,1.15,0,0.1,-0.1,-1.2];
        Source.P = 100*ones(size(Source.R));  % Watts

    case 2
        pfc.R = [-1, 0,   0,  -1, - 1,0.5,0.5, 1,1,0.1,0.1, -0.1,-0.1,-1];
        pfc.Z = [-1,-1,-1.1,-1.1,-1.2,-1.2,-1,-1,1,1,1.5,1.5,1, 1];
% Simple routine to cleanup PFC (checks for duplicate points and closes
% contour)
pfc = cleanupPFC(pfc);
        Source.R = [-0.5,0,0.5];
        Source.Z = [0,1.2,0];
        Source.P = [100,100,100];  % Watts
    case 3
        % sources from SOLPS

        % fName = 'C:\Users\jjl\Dropbox (ORNL)\SPARC\xportScan\v2y_D+Ne_INFUSE\baserun\vvfile_v2y.ogr';
        % data = dlmread(fName);
        % pfc.R = data(:,1)./1e3;
        % pfc.Z = data(:,2)./1e3;


        %%
        run_path = 'C:\Users\jjl\Dropbox (ORNL)\SPARC\xportScan\v2y_D+Ne_INFUSE\P10MW_nBC_12e19_ss_chix4';
        Case = load_solps_case(run_path,1,[],'prad');







        %%
        % figure; hold on; box on; grid on; set(gcf,'color','w'); set(gca,'fontsize',14,'fontweight','bold');
        % plot(data(:,1)./1e3,data(:,2)./1e3,'o-')
        % mesh = plot_mesh_extra(run_path,0,1);
        %
        % for i = 1:size(mesh.r,2)-1
        %     rr(i) = mesh.r(1,i);
        %     zz(i) = mesh.z(1,i);
        % end
        %
        % plot(rr,zz,'rx-')


        nUse = numel(Case.Geo.r2d_cen);
        % nUse = 220;
        Source.R = Case.Geo.r2d_cen(1:nUse);
        Source.Z = Case.Geo.z2d_cen(1:nUse);
        Source.P = Case.Prad.pRad2D(1:nUse);

        

end

%% build initial PFC from mesh.extra
 mesh = plot_mesh_extra(run_path);
[pfc.R,pfc.Z] = cleanup_mesh_extra(mesh);

% Simple routine to cleanup PFC (checks for duplicate points and closes
% contour)
pfc = cleanupPFC(pfc,0);





%% Bring in SOLPS grid points

for k = [1,2]
    if k == 1
        rPtsAdd = Case.Geo.crx(2,2:end,1);
        zPtsAdd = Case.Geo.cry(2,2:end,1);
    else
        rPtsAdd = Case.Geo.crx(end,2:end,1);
        zPtsAdd = Case.Geo.cry(end,2:end,1);
    end

    minD = nan(size(rPtsAdd));
    indD = nan(size(rPtsAdd));
    for j = 1:length(rPtsAdd)
        p0 = [rPtsAdd(j),zPtsAdd(j)];

        for i = 1:pfc.nSeg
            p1 = [pfc.R(i),pfc.Z(i)];
            p2 = [pfc.R(i+1),pfc.Z(i+1)];
            [dist(i),pu,CUTOFF,pu_unchecked,dist_unchecked] = distance_point_to_line_seg(p1,p2,p0);
        end
        [minD(j),indD(j)] = min(dist);

        % i = indD;
        % p1 = [pfc.R(i),pfc.Z(i)];
        % p2 = [pfc.R(i+1),pfc.Z(i+1)];
        % plot([p1(1),p2(1)],[p1(2),p2(2)],'linew',4)

    end

    indsToRemove = unique(indD);

    if max(minD) > 1e-3
        error('Grid points do not match PFC well enough')
    end

    % Add points into new array
    if k == 1
        Rnew = pfc.R(1:min(indD));
        Znew = pfc.Z(1:min(indD));
        indD1 = indD;
    else
        Rnew = [Rnew,pfc.R(max(indD1)+1:min(indD))];
        Znew = [Znew,pfc.Z(max(indD1)+1:min(indD))];
    end

    % find which B2 point is closest
    [~,indClose] = min(sqrt((Rnew(end) - rPtsAdd).^2 + (Znew(end) - zPtsAdd).^2));

    if indClose ~= 1
        rPtsAdd = fliplr(rPtsAdd);
        zPtsAdd = fliplr(zPtsAdd);
    end
    
    indRef(k) = length(Rnew) + 1;
    Rnew = [Rnew,rPtsAdd];
    Znew = [Znew,zPtsAdd];
end

Rnew = [Rnew,pfc.R(max(indD)+1:end)];
Znew = [Znew,pfc.Z(max(indD)+1:end)];

figure; hold on; box on; grid on; set(gcf,'color','w'); set(gca,'fontsize',14,'fontweight','bold');
plot(pfc.R,pfc.Z,'k-')
plot(Rnew,Znew,'r.-')

pfc.R = Rnew;
pfc.Z = Znew;

%% Add segment length and area
pfc.nSeg = length(pfc.R) - 1;
pfc.L = sqrt(diff(pfc.R).^2 + diff(pfc.Z).^2);

% calculate area
for i = 1:pfc.nSeg
    R1 = pfc.R(i);    R2 = pfc.R(i+1);
    Z1 = pfc.Z(i);    Z2 = pfc.Z(i+1);
    pfc.A(i) = 2*pi*sqrt( (R2 - R1)^2 + (Z2 - Z1)^2 )*mean([R1,R2]);
end


%% if source is outside pfc toss it
pTol = 1e-3;

nOrig = length(Source.P);
iKeep = Source.P >= pTol;
Source.R = Source.R(iKeep);
Source.Z = Source.Z(iKeep);
Source.P = Source.P(iKeep);
fprintf('Dropped %d source points out of %d from power threshold\n',nOrig-length(Source.P),nOrig)

iKeep = inpolygon(Source.R,Source.Z,pfc.R,pfc.Z);
nOrig = length(Source.P);
Source.R = Source.R(iKeep);
Source.Z = Source.Z(iKeep);
Source.P = Source.P(iKeep);
fprintf('Dropped %d source points out of %d from inside threshold\n',nOrig-length(Source.P),nOrig)



%%
%--------------------------------------------------------------------------
% Define rays.
% A series of rays are checked isotropically from each source point
nRays = 200;   % Number of rays per point
rayLength = 3;  % Length of each ray (must result in an intersection)
% rayPts = 1000;  % Resolution along each ray to check for intersections. Too coarse may miss grazing rays.

dTheta = 2*pi/nRays;
theta = 0:dTheta:2*pi - dTheta;



%%
cSeries = turbo;
figure; hold on; box on; set(gcf,'color','w');set(gca,'fontsize',14,'fontweight','bold'); grid on; axis equal;
plot(pfc.R,pfc.Z,'k-')
% plot(Source.R,Source.Z,'ko')

vMin = 0;
vMax = max(Source.P);

indColor = round((Source.P-vMin)./(vMax-vMin)*255 + 1);
indColor = max(indColor,1);
indColor = min(indColor,256);

scatter(Source.R(Source.P > 1e-6),Source.Z(Source.P > 1e-6),32,cSeries(indColor((Source.P > 1e-6)),:),'filled','marker','o')
% scatter(x,y,42,cSeries(i,:),'filled','marker','o')

%%
DEBUG = 0;

count = zeros(1,pfc.nSeg);
power = zeros(1,pfc.nSeg);
for iPt = 1:length(Source.R)

    if mod(iPt,100) == 0
    fprintf('Working on point %d of %d\n',iPt,length(Source.R))
    end
    rStart = Source.R(iPt);
    zStart = Source.Z(iPt);


    for iRay = 1:nRays
        rEnd = rStart + rayLength*cos(theta(iRay));
        zEnd = zStart + rayLength*sin(theta(iRay));

        if METHOD == 1
            % Define ray
            rRay = linspace(rStart,rEnd,rayPts);
            zRay = linspace(zStart,zEnd,rayPts);

            % Check for intersections with pfcs
            in = inpolygon(rRay,zRay,pfc.R,pfc.Z);
            ind = find(~in,1,'first') - 1;

            % If no intersections ray length must be increased
            if isempty(ind)
                plot(rRay,zRay)
                error('Increase ray length or something')
            end

            % Find precise intersection point
            p1 = [rRay(ind),zRay(ind)];
            p2 = [rRay(ind+1),zRay(ind+1)];
            [pint1,ierr,found_ind,int_count] = int_line_curve(p1,p2,pfc.R,pfc.Z);

            % Can be a problem if pfc is overlapping or ray is too coarse
            if int_count ~= 1
                for j = 1:int_count
                    plot(pint1(j,1),pint1(j,2),'x')
                end
                error('More than one intersection? Too detailed PFC or too coarse ray? Or PFC structure is not neatly ordered')
            end
        elseif METHOD == 2

            [pint1,ierr,found_ind,int_count] = int_line_curve([rStart,zStart],[rEnd,zEnd],pfc.R,pfc.Z);
            Ltest = sqrt( (rStart - pint1(:,1)).^2 + (zStart - pint1(:,2)).^2 );
            [~,iUse] = min(Ltest);
            pint1 = pint1(iUse,:);
            found_ind = found_ind(iUse);
        end


        % Tally counts and power
        count(found_ind) = count(found_ind) + 1;
        power(found_ind) = power(found_ind) + Source.P(iPt)/nRays;

        if DEBUG
            % plot(rRay,zRay)
            % plot(rRay(ind),zRay(ind),'x')
            plot(pint1(1),pint1(2),'x')
        end
    end
end


toc

%%  Make a plot
c = parula;
n = size(c,1);

if PLOT_TYPE == 1
    w = power; % Total power per element
    title('Power (W) per element')
elseif PLOT_TYPE == 2
    w = power./pfc.L; % Power density per element
    title('Linear power density (W/m) per element')
elseif PLOT_TYPE == 3
    w = power./pfc.A; % Power density per element
    title('Areal power density (W/m^2) per element')
end

% These are color bar ranges
% cMax = Source.P;
cMax = max(w);
cMin = 0;

p = (w - cMin)./(cMax - cMin); % Relative position in color table
cInd = floor(p*(n-1)) + 1;
cInd(cInd<1) = 1;
cInd(cInd>n) = n;
cInd(isnan(cInd)) = 1; % Sometimes L = 0

for i = 1:pfc.nSeg
    plot(pfc.R(i:i+1),pfc.Z(i:i+1),'-','linew',2,'color',c(cInd(i),:))
end
colorbar;
clim([cMin,cMax])

plot(pfc.R,pfc.Z,'k.','markersize',8)


figure; hold on; box on; grid on; set(gcf,'color','w'); set(gca,'fontsize',14,'fontweight','bold');
plot(Case.DsData.dsr(2:end-1),w(indRef(2):indRef(2)+35))

%% Cleanup routine
function pfc = cleanupPFC(pfc,closeIt)
if nargin < 2
    closeIt = 1;
end

if closeIt
    %% Add first point to make closed contour
    if (pfc.R(1) == pfc.R(end)) && (pfc.Z(1) == pfc.Z(end))
    else
        pfc.R(end+1) = pfc.R(1);
        pfc.Z(end+1) = pfc.Z(1);
    end
end

pfc.nSeg = length(pfc.R) - 1;

%% remove repeated points (only works for one set of adjacent points)
iCount = 1;
for i = 1:pfc.nSeg
    if (pfc.R(i) == pfc.R(i+1)) && (pfc.Z(i) == pfc.Z(i+1))
        continue
    end
    rNew(iCount) = pfc.R(i);
    zNew(iCount) = pfc.Z(i);
    iCount = iCount + 1;
end
rNew(end+1) = pfc.R(end);
zNew(end+1) = pfc.Z(end);
pfc.R = rNew;
pfc.Z = zNew;
end

