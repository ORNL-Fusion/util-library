clearvars;


config_name = {'D3-6'};

verbose = 1;
SIMPLIFY_COILS = 1;
PLOT_COILS = 1;
PLOT_LC = 'right';
PLOT_SCATTER = 0;
PLOT_TRIANGULATION = 1;
TRIANGULATION_ALPHA = []; % 0.05, 0.03, 0.08 --  bigger smoother - smaller tighter


Geo = get_MPEX_geometry;

%%
for i = 1:length(config_name)
    figure; set(gcf,'color','w'); box on; grid on; hold on; set(gca,'fontsize',14)

    %% Setup coils and bfield for this config
    if SIMPLIFY_COILS
        [Coil,windingCurrent,CoilGeometry,currentPerWinding] = build_MPEX_coils_jackson_hybrid(config_name{i},[],3,3,verbose);
    else
        [Coil,windingCurrent,CoilGeometry,currentPerWinding] = build_MPEX_coils_jackson(config_name{i},verbose);
    end
    Bfield.coil = Coil; Bfield.current = windingCurrent; Bfield.type = 'MPEX';


    %% Plot coils and vessel
    if PLOT_COILS
        [rcoil,zcoil] = get_coil_cross_sections(CoilGeometry);
        plot_coil_cross_section(rcoil,zcoil,0,currentPerWinding);
    else
        rcoil = [];
        zcoil = [];
    end
    plot(Geo.Vessel.z,Geo.Vessel.r,'k-','LineWidth',2)
    plot(Geo.Target.z,Geo.Target.r,'k-','LineWidth',5)
    xlabel('Z (m)')
    ylabel('R (m)')
    axis([-4,9,0,0.5])

    %% Find LUFS
    need_lufs_overlay = PLOT_SCATTER || PLOT_TRIANGULATION;
    if need_lufs_overlay
        f_lufs = find_lufs_MPEX(Bfield,Geo);
    else
        f_lufs = struct('r',[],'z',[],'hit_rz',[]);
    end
    

    %% Now let's compute connection length in space
    % nlines_scattered = 50;
    % Rmin_scattered = 0.3;    
    % Rmax_scattered = 1;
    % RSeed = linspace(Rmin_scattered,Rmax_scattered,nlines_scattered);
    RSeed = unique([linspace(0,0.1,50),linspace(0.1,0.3,150),linspace(0.3,1,150)]);
    Zmin = min(Geo.Vessel.z);
    Zmax = max(Geo.Vessel.z);
    dz_fieldline0 = 0.01;


    need_left = strcmp(PLOT_LC,'left') || strcmp(PLOT_LC,'total');
    need_right = strcmp(PLOT_LC,'right') || strcmp(PLOT_LC,'total');
    need_total = strcmp(PLOT_LC,'total');
    
    if strcmp(PLOT_LC,'left')
        plot_title = 'Left connection length (m)';
    elseif strcmp(PLOT_LC,'right')
        plot_title = 'Right connection length (m)';
    elseif strcmp(PLOT_LC,'total')
        plot_title = 'Total connection length (m)';
    else
        error('Unknown PLOT_LC %s',PLOT_LC)
    end
    
    [r_plot,z_plot,L_forward_plot,L_backward_plot,L_total_plot] = calc_connection_length_scattered(Bfield,Geo,RSeed,Zmin,Zmax,dz_fieldline0,need_left,need_right,need_total);
    
    if strcmp(PLOT_LC,'left')
        L_plot = L_backward_plot;
    elseif strcmp(PLOT_LC,'right')
        L_plot = L_forward_plot;
    else
        L_plot = L_total_plot;
    end
    
    if PLOT_SCATTER
        figure; set(gcf,'color','w')
        scatter(z_plot,r_plot,10,L_plot,'filled')
        hold on
        plot(Geo.Vessel.z,Geo.Vessel.r,'k-','LineWidth',2)
        plot(Geo.Target.z,Geo.Target.r,'k-','LineWidth',5)
        box on; grid on; set(gca,'fontsize',14)
        xlabel('Z (m)')
        ylabel('R (m)')
        title(plot_title)
        colorbar
        axis([-4,9,0,0.5])
        if PLOT_COILS
            plot_coil_cross_section(rcoil,zcoil,0);
        end
        if ~isempty(f_lufs.z)
            plot(f_lufs.z,f_lufs.r,'r')
        end
    end

    if PLOT_TRIANGULATION
        figure; set(gcf,'color','w')
        plot_triangulated_map(z_plot,r_plot,L_plot,Geo,plot_title,TRIANGULATION_ALPHA)
        if PLOT_COILS
            plot_coil_cross_section(rcoil,zcoil,0);
        end
        if ~isempty(f_lufs.z)
            plot(f_lufs.z,f_lufs.r,'r')
        end
    end

end

function f = find_lufs_MPEX(Bfield,Geo,RMax)

% if nargin < 3
    num_lines = 10;
% end
% if nargin < 5
    nbisect = 8;
% end

ZMin = min(Geo.Vessel.z) + 0.01;  % Below this Z lines are not checked for cutoff
if nargin < 3
    RMax = max(Geo.Target.r) + 0.01;
end
RMin = 1e-3;

% set up fieldline following
L = Geo.Target.z(1) - ZMin;
dz_fieldline = -0.01; 
nsteps = round(abs(L/dz_fieldline)); 
dz_fieldline = sign(dz_fieldline)*L/nsteps;

rr = linspace(RMin,RMax,num_lines);
phistart = zeros(size(rr));
fl = follow_fieldlines_rzphi_dz(Bfield,rr,Geo.Target.z(1),phistart,dz_fieldline,nsteps);
fl = clip_fl_at_vessel(fl,Geo.Vessel.r,Geo.Vessel.z);
isgood = ~isnan(fl.z(end,:));
ind_last_good_line = find(isgood,1,'last');

if isempty(ind_last_good_line)
    error('Could not identify LUFS bracket: all coarse lines hit the vessel before reaching ZMin.')
end
if ind_last_good_line == numel(rr)
    error('Could not identify LUFS bracket: all coarse lines passed through. Increase RMax.')
end

r_good = rr(ind_last_good_line);
r_bad = rr(ind_last_good_line+1);
fl_good.r = fl.r(:,ind_last_good_line);
fl_good.z = fl.z(:,ind_last_good_line);
fl_good.phi = fl.phi(:,ind_last_good_line);
fl_bad.r = fl.r(:,ind_last_good_line+1);
fl_bad.z = fl.z(:,ind_last_good_line+1);
fl_bad.phi = fl.phi(:,ind_last_good_line+1);

for ib = 1:nbisect
    r_mid = 0.5*(r_good + r_bad);
    fl_mid = trace_lufs_line(Bfield,Geo,r_mid,Geo.Target.z(1),dz_fieldline,nsteps);
    if ~isnan(fl_mid.z(end))
        r_good = r_mid;
        fl_good = fl_mid;
    else
        r_bad = r_mid;
        fl_bad = fl_mid;
    end
end

ind_hit = find(~isnan(fl_bad.z),1,'last');
if isempty(ind_hit)
    hit_rz = [NaN,NaN];
else
    hit_rz = [fl_bad.r(ind_hit),fl_bad.z(ind_hit)];
end

f.r = fl_good.r;
f.z = fl_good.z;
f.hit_rz = hit_rz;


end


function fl_single = trace_lufs_line(Bfield,Geo,r0,z0,dz_fieldline,nsteps)
fl_raw = follow_fieldlines_rzphi_dz(Bfield,r0,z0,0,dz_fieldline,nsteps);
fl_single = clip_fl_at_vessel(fl_raw,Geo.Vessel.r,Geo.Vessel.z);
end



function fout = clip_fl_at_vessel(f,vessel_r,vessel_z)
% Clip field lines at the first point that leaves the vessel
% Same field line structure returned, but points after vessel exit are set
% to NaN
fout.r = NaN(size(f.r));
fout.z = NaN(size(f.z));
fout.phi = NaN(size(f.phi));

for i = 1:size(f.r,2)
    isin2 = inpolygon(f.r(:,i),f.z(:,i),vessel_r,vessel_z);
    % Keep points up to the last point still inside the vessel
    is2 = find(isin2 == 0,1,'first') - 1;
    if isempty(is2)
        is2 = length(isin2);
    end
    if is2 > 0
        fout.r(1:is2,i) = f.r(1:is2,i);
        fout.z(1:is2,i) = f.z(1:is2,i);
        fout.phi(1:is2,i) = f.phi(1:is2,i);
    end
end
end


function fout = mask_fl_at_vessel(f,vessel_r,vessel_z)
% Same field line structure returned, but points outside the vessel are NaN
fout = f;
isin = inpolygon(f.r,f.z,vessel_r,vessel_z);
fout.r(~isin) = NaN;
fout.z(~isin) = NaN;
if isfield(f,'phi')
    fout.phi(~isin) = NaN;
end
end


function [Lc_forward,Lc_backward,Lc_total] = calc_connection_length_sections(f,calc_forward,calc_backward,calc_total)
% Compute connection length within each continuous in-vessel section
if nargin < 2
    calc_forward = 1;
end
if nargin < 3
    calc_backward = 1;
end
if nargin < 4
    calc_total = 1;
end

Lc_forward = NaN(size(f.r));
Lc_backward = NaN(size(f.r));
Lc_total = NaN(size(f.r));

for i = 1:size(f.r,2)
    isgood = ~isnan(f.r(:,i)) & ~isnan(f.z(:,i));
    if ~any(isgood)
        continue
    end

    idiff = diff([0;isgood;0]);
    istart = find(idiff == 1);
    iend = find(idiff == -1) - 1;
    for j = 1:length(istart)
        ind = istart(j):iend(j);
        Lseg = [0;cumsum(sqrt(diff(f.r(ind,i)).^2 + diff(f.z(ind,i)).^2))];
        if calc_backward
            Lc_backward(ind,i) = Lseg;
        end
        if calc_forward
            Lc_forward(ind,i) = Lseg(end) - Lseg;
        end
        if calc_total
            Lc_total(ind,i) = Lseg(end);
        end
    end
end
end


function [r_plot,z_plot,L_forward_plot,L_backward_plot,L_total_plot] = calc_connection_length_scattered(Bfield,Geo,RSeed,Zmin,Zmax,dz_fieldline0,need_left,need_right,need_total)
% RSeed = linspace(Rmin_scattered,Rmax_scattered,nlines_scattered);
ZSeed = Zmin*ones(size(RSeed));
L = Zmax - Zmin;
nsteps = max(1,ceil(abs(L/dz_fieldline0)));
dz_fieldline = L/nsteps;
phistart = zeros(size(RSeed));
fl_all = follow_fieldlines_rzphi_dz(Bfield,RSeed,ZSeed(1),phistart,dz_fieldline,nsteps);
fl_all = mask_fl_at_vessel(fl_all,Geo.Vessel.r,Geo.Vessel.z);
[Lc_forward_all,Lc_backward_all,Lc_total_all] = calc_connection_length_sections(fl_all,need_right,need_left,need_total);

r_plot = [];
z_plot = [];
L_forward_plot = [];
L_backward_plot = [];
L_total_plot = [];
for iline = 1:size(fl_all.r,2)
    if need_left && ~need_right && ~need_total
        ivalid = find(~isnan(Lc_backward_all(:,iline)));
    elseif need_right && ~need_left && ~need_total
        ivalid = find(~isnan(Lc_forward_all(:,iline)));
    else
        ivalid = find(~isnan(Lc_total_all(:,iline)));
    end
    if isempty(ivalid)
        continue
    end
    rline = fl_all.r(ivalid,iline);
    zline = fl_all.z(ivalid,iline);

    r_plot = [r_plot;rline(:)];
    z_plot = [z_plot;zline(:)];
    L_forward_plot = [L_forward_plot;Lc_forward_all(ivalid,iline)];
    L_backward_plot = [L_backward_plot;Lc_backward_all(ivalid,iline)];
    L_total_plot = [L_total_plot;Lc_total_all(ivalid,iline)];
end
end


function plot_triangulated_map(z_plot,r_plot,L_plot,Geo,plot_title_str,triangulation_alpha)
dt = delaunayTriangulation(z_plot(:),r_plot(:));
tri = dt.ConnectivityList;
hold on
ztri = reshape(dt.Points(tri,1),size(tri));
rtri = reshape(dt.Points(tri,2),size(tri));
ztri_centroid = mean(ztri,2);
rtri_centroid = mean(rtri,2);

if isempty(triangulation_alpha)
    shp = alphaShape(z_plot(:),r_plot(:));
    shp.Alpha = criticalAlpha(shp,'one-region');
else
    shp = alphaShape(z_plot(:),r_plot(:),triangulation_alpha);
end

inside_alpha = inShape(shp,ztri_centroid,rtri_centroid);
inside_vessel = inpolygon(ztri(:,1),rtri(:,1),Geo.Vessel.z,Geo.Vessel.r) & ...
    inpolygon(ztri(:,2),rtri(:,2),Geo.Vessel.z,Geo.Vessel.r) & ...
    inpolygon(ztri(:,3),rtri(:,3),Geo.Vessel.z,Geo.Vessel.r) & ...
    inpolygon(ztri_centroid,rtri_centroid,Geo.Vessel.z,Geo.Vessel.r);
inside_tri = inside_alpha & inside_vessel;
tri = tri(inside_tri,:);

patch('Faces',tri, ...
    'Vertices',dt.Points, ...
    'FaceVertexCData',L_plot(:), ...
    'FaceColor','interp', ...
    'EdgeColor','none');
plot(Geo.Vessel.z,Geo.Vessel.r,'k-','LineWidth',2)
plot(Geo.Target.z,Geo.Target.r,'k-','LineWidth',5)
box on; grid on; set(gca,'fontsize',14)
xlabel('Z (m)')
ylabel('R (m)')
title(plot_title_str)
colorbar
axis([-4,9,0,0.5])
end

