clearvars;

PLOT_TYPE = 3; % 1 is power, 2 is power density (W/m), 3 power density (W/m^2)

%%
%--------------------------------------------------------------------------
% Define PFC and Source
% Note: PFC must be non-overlapping and simply sorted 
%--------------------------------------------------------------------------

CASE = 1; % 1 = tokamak structure from gfile, 2 = simple defined geometry

switch CASE
    case 1
        gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Haskey\g179477.02280';
        g = readg_g3d(gfile_name);
        % remove weird point specific to this gfile
        g.lim(:,7) = [];
        pfc.R = g.lim(1,:);
        pfc.Z = g.lim(2,:);
        
        Source.R = [1.5,1.45,2.3,2.3,2.3,1.5];
        Source.Z = [0,1.15,0,0.1,-0.1,-1.2];
        Source.P = 100*ones(size(Source.R));  % Watts

    case 2
        pfc.R = [-1, 0,   0,  -1, - 1,0.5,0.5, 1,1,0.1,0.1, -0.1,-0.1,-1];
        pfc.Z = [-1,-1,-1.1,-1.1,-1.2,-1.2,-1,-1,1,1,1.5,1.5,1, 1];

        Source.R = [-0.5,0,0.5];
        Source.Z = [0,1.2,0];
        Source.P = [100,100,100];  % Watts
        
end

% Simple routine to cleanup PFC (checks for duplicate points and closes
% contour)
pfc = cleanupPFC(pfc);
% Add segment length and area
pfc.nSeg = length(pfc.R) - 1;
pfc.L = sqrt(diff(pfc.R).^2 + diff(pfc.Z).^2);

% calculate area 
for i = 1:pfc.nSeg
    R1 = pfc.R(i);    R2 = pfc.R(i+1);
    Z1 = pfc.Z(i);    Z2 = pfc.Z(i+1);    
    pfc.A(i) = 2*pi*sqrt( (R2 - R1)^2 + (Z2 - Z1)^2 )*mean([R1,R2]);
end

%%
%--------------------------------------------------------------------------
% Define rays. 
% A series of rays are checked isotropically from each source point
nRays = 1000;   % Number of rays per point
rayLength = 3;  % Length of each ray (must result in an intersection)
rayPts = 1000;  % Resolution along each ray to check for intersections. Too coarse may miss grazing rays.

dTheta = 2*pi/nRays;
theta = 0:dTheta:2*pi - dTheta;



%%
figure; hold on; box on; set(gcf,'color','w');set(gca,'fontsize',14,'fontweight','bold'); grid on; axis equal;
plot(pfc.R,pfc.Z,'k-')
plot(Source.R,Source.Z,'ko')

DEBUG = 0;

count = zeros(1,pfc.nSeg);
power = zeros(1,pfc.nSeg);
for iPt = 1:length(Source.R)
    rStart = Source.R(iPt);
    zStart = Source.Z(iPt);

    for iRay = 1:nRays
        rEnd = rStart + rayLength*cos(theta(iRay));
        zEnd = zStart + rayLength*sin(theta(iRay));

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

        % Tally counts and power
        count(found_ind) = count(found_ind) + 1;
        power(found_ind) = power(found_ind) + Source.P(iPt)/nRays;

        if DEBUG
            plot(rRay,zRay)
            plot(rRay(ind),zRay(ind),'x')
            plot(pint1(1),pint1(2),'x')
        end
    end
end


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


%% Cleanup routine
function pfc = cleanupPFC(pfc)
%% Add first point to make closed contour
if (pfc.R(1) == pfc.R(end)) && (pfc.Z(1) == pfc.Z(end))
else
    pfc.R(end+1) = pfc.R(1);
    pfc.Z(end+1) = pfc.Z(1);
end
pfc.nSeg = length(pfc.R) - 1;

%%
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

