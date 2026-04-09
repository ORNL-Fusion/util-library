clearvars;


config_name = {'D3-6'};

verbose = 1;
SIMPLIFY_COILS = 1;
PLOT_COILS = 1;
% METHOD = 'structured';
METHOD = 'scattered';
PLOT_LC = 'right';


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
    end
    plot(Geo.Vessel.z,Geo.Vessel.r,'k-','LineWidth',2)
    plot(Geo.Target.z,Geo.Target.r,'r-','LineWidth',2)
    xlabel('Z (m)')
    ylabel('R (m)')
    axis([-4,9,0,0.5])

    %% Find LUFS
    % f_lufs = find_lufs_MPEX(Bfield,Geo);
    % plot(f_lufs.z,f_lufs.r)

    %% Now let's compute connection length in space
    nz = 1000;
    nr = 100;
    nlines_scattered = 200;
    Rmin_structured = 0.2;
    Rmax_structured = 0.3;
    Rmin_scattered = 0;
    Rmax_scattered = 0.7;
    Zmin = min(Geo.Vessel.z);
    Zmax = max(Geo.Vessel.z);
    dz_fieldline0 = 0.01;
    dmax_interp = 0.10;
    hull_shrink_factor = 0.2;
    need_left = strcmp(PLOT_LC,'left') || strcmp(PLOT_LC,'total');
    need_right = strcmp(PLOT_LC,'right') || strcmp(PLOT_LC,'total');
    
    if strcmp(PLOT_LC,'left')
        plot_title = 'Left connection length';
    elseif strcmp(PLOT_LC,'right')
        plot_title = 'Right connection length';
    elseif strcmp(PLOT_LC,'total')
        plot_title = 'Total connection length';
    else
        error('Unknown PLOT_LC %s',PLOT_LC)
    end
    
    if strcmp(METHOD,'structured')
        Redges = linspace(Rmin_structured,Rmax_structured,nr+1);
        Zedges = linspace(Zmin,Zmax,nz+1);
        R1d = 0.5*(Redges(1:end-1) + Redges(2:end));
        Z1d = 0.5*(Zedges(1:end-1) + Zedges(2:end));
        [R2d,Z2d] = ndgrid(R1d,Z1d);
        [Lc2d_forward,Lc2d_backward,Lc2d_total] = calc_connection_length_structured(Bfield,Geo,R2d,Z2d,Zmin,dz_fieldline0,need_left,need_right,strcmp(PLOT_LC,'total'));
        
        if strcmp(PLOT_LC,'left')
            Lc2d_plot = Lc2d_backward;
        elseif strcmp(PLOT_LC,'right')
            Lc2d_plot = Lc2d_forward;
        else
            Lc2d_plot = Lc2d_total;
        end
        
        figure; set(gcf,'color','w')
        plot_cell_centered_patches(Z2d,R2d,Lc2d_plot,Geo,['Structured ',plot_title])
    elseif strcmp(METHOD,'scattered')
        [R2d,Z2d,Lc2d_forward,Lc2d_backward,Lc2d_total,r_plot,z_plot,L_forward_plot,L_backward_plot,L_total_plot] = calc_connection_length_scattered(Bfield,Geo,Rmin_scattered,Rmax_scattered,Zmin,Zmax,nlines_scattered,nr,nz,dz_fieldline0,dmax_interp,hull_shrink_factor,need_left,need_right,strcmp(PLOT_LC,'total'));
        
        if strcmp(PLOT_LC,'left')
            L_plot = L_backward_plot;
            Lc2d_plot = Lc2d_backward;
        elseif strcmp(PLOT_LC,'right')
            L_plot = L_forward_plot;
            Lc2d_plot = Lc2d_forward;
        else
            L_plot = L_total_plot;
            Lc2d_plot = Lc2d_total;
        end
        
        figure; set(gcf,'color','w')
        scatter(z_plot,r_plot,10,L_plot,'filled')
        hold on
        plot(Geo.Vessel.z,Geo.Vessel.r,'k-','LineWidth',2)
        plot(Geo.Target.z,Geo.Target.r,'r-','LineWidth',2)
        box on; grid on; set(gca,'fontsize',14)
        xlabel('Z (m)')
        ylabel('R (m)')
        title([plot_title,' from dense field line samples'])
        colorbar
        axis([-4,9,0,0.5])
        
        figure; set(gcf,'color','w')
        plot_cell_centered_patches(Z2d,R2d,Lc2d_plot,Geo,['Interpolated ',plot_title])
    else
        error('Unknown METHOD %s',METHOD)
    end


end


%%
function f = find_lufs_MPEX(Bfield,Geo,RMax)

% if nargin < 3
    num_lines = 10;
% end
% if nargin < 5
    nresolve = 3;
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

for iresolve = 1:nresolve
    
    % Forward and reverse lines from target
    if iresolve == 1
        rr = linspace(RMin,RMax,num_lines);
    else
        rr = linspace(rr(ind_last_good_line-1),rr(ind_last_good_line+1),num_lines*iresolve);
    end
    zz = Geo.Target.z(1)*ones(size(rr));

    phistart = zeros(size(rr));    
    fl = follow_fieldlines_rzphi_dz(Bfield,rr,zz(1),phistart,dz_fieldline,nsteps);
    fl = clip_fl_at_vessel(fl,Geo.Vessel.r,Geo.Vessel.z);
    ind_last_good_line = find(~isnan(fl.z(end,:)),1,'last');
    if iresolve == nresolve
        ind_hit = find(~isnan(fl.z(:,ind_last_good_line+1)),1,'last');
        zhit = fl.z(ind_hit,ind_last_good_line+1);
        rhit = fl.r(ind_hit,ind_last_good_line+1);
        hit_rz = [rhit,zhit];
    end
    
    if isempty(ind_last_good_line)        
        plotit = 1;
    else
        plotit = 0;
    end
    
    if plotit
        if iresolve == 1
            figure; hold on; box on;
            plot(fl.z,fl.r,'linewidth',2);
            set(gca,'fontsize',14)
            xlabel('Z [m]','fontsize',14)
            ylabel('R [m]','fontsize',14)
            plot(Geo.Vessel.z,Geo.Vessel.r,'k-','LineWidth',2)
            plot(Geo.Target.z,Geo.Target.r,'r-','LineWidth',2)

        else
            plot(fl.z,fl.r,'linewidth',2,'color',cf(iresolve,:));
            plot(fl.z(:,ind_last_good_line),fl.r(:,ind_last_good_line),'r')
        end
        if isempty(ind_last_good_line)
            error('Could not identify lcfs: all lines passed through')
        end
    end
end
f.r = fl.r(:,ind_last_good_line);
f.z = fl.z(:,ind_last_good_line);
f.hit_rz = hit_rz;


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


function [Lc2d_forward,Lc2d_backward,Lc2d_total] = calc_connection_length_structured(Bfield,Geo,R2d,Z2d,Zmin,dz_fieldline0,need_left,need_right,need_total)
Lc2d_forward = NaN(size(R2d));
Lc2d_backward = NaN(size(R2d));
Lc2d_total = NaN(size(R2d));

for ir = 1:size(R2d,1)
    fprintf('Working on ir %d of %d\n',ir,size(R2d,1))
    for jz = 1:size(R2d,2)
        if need_right
            L = Geo.Target.z(1) - Z2d(ir,jz);
            nsteps = max(1,ceil(abs(L/dz_fieldline0)));
            dz_fieldline = L/nsteps;
            fl = follow_fieldlines_rzphi_dz(Bfield,R2d(ir,jz),Z2d(ir,jz),0,dz_fieldline,nsteps);
            fl_mask = mask_fl_at_vessel(fl,Geo.Vessel.r,Geo.Vessel.z);
            [Lc_forward_fl,~,~] = calc_connection_length_sections(fl_mask,1,0,0);
            Lc2d_forward(ir,jz) = Lc_forward_fl(1);
        end

        if need_left
            L = Zmin - Z2d(ir,jz);
            nsteps = max(1,ceil(abs(L/dz_fieldline0)));
            dz_fieldline = L/nsteps;
            fl = follow_fieldlines_rzphi_dz(Bfield,R2d(ir,jz),Z2d(ir,jz),0,dz_fieldline,nsteps);
            fl_mask = mask_fl_at_vessel(fl,Geo.Vessel.r,Geo.Vessel.z);
            [~,Lc_backward_fl,~] = calc_connection_length_sections(fl_mask,0,1,0);
            Lc2d_backward(ir,jz) = Lc_backward_fl(1);
        end

        if need_total
            Lc2d_total(ir,jz) = Lc2d_forward(ir,jz) + Lc2d_backward(ir,jz);
        end
    end
end
end


function [R2d,Z2d,Lc2d_forward,Lc2d_backward,Lc2d_total,r_plot,z_plot,L_forward_plot,L_backward_plot,L_total_plot] = calc_connection_length_scattered(Bfield,Geo,Rmin_scattered,Rmax_scattered,Zmin,Zmax,nlines_scattered,nr,nz,dz_fieldline0,dmax_interp,hull_shrink_factor,need_left,need_right,need_total)
RSeed = linspace(Rmin_scattered,Rmax_scattered,nlines_scattered);
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

Redges = linspace(Rmin_scattered,Rmax_scattered,nr+1);
Zedges = linspace(Zmin,Zmax,nz+1);
R1d = 0.5*(Redges(1:end-1) + Redges(2:end));
Z1d = 0.5*(Zedges(1:end-1) + Zedges(2:end));
[R2d,Z2d] = ndgrid(R1d,Z1d);
outside_vessel = ~inpolygon(R2d,Z2d,Geo.Vessel.r,Geo.Vessel.z);
ihull = boundary(z_plot,r_plot,hull_shrink_factor);
inside_data_hull = inpolygon(Z2d,R2d,z_plot(ihull),r_plot(ihull));
Dmin = NaN(size(R2d));
for ir = 1:size(R2d,1)
    for jz = 1:size(R2d,2)
        if ~outside_vessel(ir,jz) && inside_data_hull(ir,jz)
            Dmin(ir,jz) = min(sqrt((r_plot - R2d(ir,jz)).^2 + (z_plot - Z2d(ir,jz)).^2));
        end
    end
end
istoo_far = Dmin > dmax_interp;
if any(istoo_far(~outside_vessel & inside_data_hull),'all')
    warning('Some in-vessel interpolation points are farther than dmax_interp = %.2f m from the field-line samples. Consider increasing nlines_scattered, increasing nr/nz for the sampled lines, or reducing the interpolation domain.',dmax_interp)
end

igood = ~outside_vessel & inside_data_hull & ~istoo_far;
Lc2d_forward = NaN(size(R2d));
Lc2d_backward = NaN(size(R2d));
Lc2d_total = NaN(size(R2d));
if need_right
    F_forward = scatteredInterpolant(r_plot,z_plot,L_forward_plot,'linear','none');
    Lc2d_forward(igood) = F_forward(R2d(igood),Z2d(igood));
end
if need_left
    F_backward = scatteredInterpolant(r_plot,z_plot,L_backward_plot,'linear','none');
    Lc2d_backward(igood) = F_backward(R2d(igood),Z2d(igood));
end
if need_total
    F_total = scatteredInterpolant(r_plot,z_plot,L_total_plot,'linear','none');
    Lc2d_total(igood) = F_total(R2d(igood),Z2d(igood));
end
end


function plot_cell_centered_patches(Z2d,R2d,C2d,Geo,plot_title_str)
z_centers = Z2d(1,:);
r_centers = R2d(:,1);
z_edges = centers_to_edges(z_centers);
r_edges = centers_to_edges(r_centers);
[Redge2d,Zedge2d] = ndgrid(r_edges,z_edges);

inside_vertices = inpolygon(Redge2d,Zedge2d,Geo.Vessel.r,Geo.Vessel.z);
inside_cells = inside_vertices(1:end-1,1:end-1) & inside_vertices(2:end,1:end-1) & ...
    inside_vertices(1:end-1,2:end) & inside_vertices(2:end,2:end);
valid_cells = isfinite(C2d) & inside_cells;

vertex_ids = reshape(1:numel(Redge2d),size(Redge2d));
v1 = vertex_ids(1:end-1,1:end-1);
v2 = vertex_ids(2:end,1:end-1);
v3 = vertex_ids(2:end,2:end);
v4 = vertex_ids(1:end-1,2:end);
faces = [v1(:), v2(:), v3(:), v4(:)];
cell_values = C2d(:);
valid_faces = valid_cells(:);

patch('Faces',faces(valid_faces,:), ...
    'Vertices',[Zedge2d(:),Redge2d(:)], ...
    'FaceVertexCData',cell_values(valid_faces), ...
    'FaceColor','flat', ...
    'EdgeColor','none');
hold on
plot(Geo.Vessel.z,Geo.Vessel.r,'k-','LineWidth',2)
plot(Geo.Target.z,Geo.Target.r,'r-','LineWidth',2)
box on; grid on; set(gca,'fontsize',14)
xlabel('Z (m)')
ylabel('R (m)')
title(plot_title_str)
colorbar
axis([-4,9,0,0.5])
end


function edges = centers_to_edges(centers)
centers = centers(:).';

if numel(centers) < 2
    error('Need at least two cell centers to infer cell edges.')
end

dcenters = diff(centers);
edges = [centers(1) - 0.5*dcenters(1), ...
    0.5*(centers(1:end-1) + centers(2:end)), ...
    centers(end) + 0.5*dcenters(end)];
end
