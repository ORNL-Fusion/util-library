clearvars;
tic;
connect_dots = 0;

gfile_name = 'C:\Work\EMC3_revival\gfiles\DIII-D\g154929.03750';
g = readg_g3d(gfile_name);

rmp.type = 'g3d';
g.toroidal_off_grid = 1;

% AXISYMMETRIC-AXISYMMETRIC-AXISYMMETRIC-AXISYMMETRIC %
% rmp.isaxisym = 1;
% rmp.current = 0;
% RMP-RMP-RMP-RMP-RMP-RMP-RMP-RMP-RMP-RMP-RMP-RMP-RMP %
% taper = [-897.705, 2362.30, 3228.76, 897.705, -2362.30, -3228.76];
% rmp = build_d3d_ccoils_jl(taper);
% rmp.type = 'g3d';
% rmp.isaxisym = 0;
% nsym = 1;
% RMP-RMP-RMP-RMP-RMP-RMP-RMP-RMP-RMP-RMP-RMP-RMP-RMP %
% taper_i = [   -3740.   3863.  -3720.   3855.  -3718.   3858. 3862.  -3791.   3884.  -3854.   3923.  -3847.];
taper_i = 3826*[-1 1 -1 1 -1 1 1 -1 1 -1 1 -1];
rmp = build_d3d_icoils_jl(taper_i); 
rmp.type = 'g3d';
rmp.isaxisym = 0;
nsym = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Rstart = 1.7;    % if [] this is set to g.rmaxis
Rend   = 2.25;    % if [] this is set to max radius of g.bdry
Rstart = 2.05;    % if [] this is set to g.rmaxis
% Rend   = 1.7;    % if [] this is set to max radius of g.bdry

Zstart = 0;
Zend = 0;
num_surfs = 10;

dphi_deg = 0.5;
num_transits = 40; % fieldlines are followed in phi each direction 

phistart_deg = 0;
%--------------------------------------------------------------------------

if isempty(Rend) 
    Rend = max(g.bdry(1,:));
end
if isempty(Rstart)
    Rstart = g.rmaxis;
end

Rstarts = linspace(Rstart,Rend,num_surfs);
Zstarts = linspace(Zstart,Zend,num_surfs);
% Zstarts = zeros(size(Rstarts));
dphi=dphi_deg*pi/180;
nsteps = floor(num_transits*2*pi/dphi);
phistart = phistart_deg*pi/180;


s1=follow_fieldlines_rzphi(g,rmp,Rstarts,Zstarts,phistart,dphi,nsteps,0);
s2=follow_fieldlines_rzphi(g,rmp,Rstarts,Zstarts,phistart,-dphi,nsteps,0);
phi = [s2.phi(end:-1:2);s1.phi];
r = [s2.r(end:-1:2,:);s1.r];
z = [s2.z(end:-1:2,:);s1.z];

% for i = 1:num_surfs
%     fprintf('Following line %i of %i\n',i,num_surfs)
%     s1=follow_fieldlines_rzphi(g,rmp,Rstarts(i),Zstarts(i),phistart,dphi,nsteps);
%     s2=follow_fieldlines_rzphi(g,rmp,Rstarts(i),Zstarts(i),phistart,-dphi,nsteps);
%     phi = [s2.phi(end:-1:2);s1.phi];
%     r(:,i) = [s2.r(end:-1:2,:);s1.r];
%     z(:,i) = [s2.z(end:-1:2,:);s1.z];
% end

if connect_dots
    for i = 1:num_surfs
        theta0 = mod(atan2(z(:,i)-g.zmaxis,r(:,i)-g.rmaxis),2*pi);    
        [~,sort_inds] = sort(theta0);
        r(:,i) = r(sort_inds,i);   
        z(:,i) = z(sort_inds,i);
        phi = phi(sort_inds);
        if ~rmp.isaxisym
            phi = mod(phi,2*pi/nsym);
            inds = find(abs(phi - phistart_deg) < 1e-6);
            rfinal(:,i) = r(inds,i);
            zfinal(:,i) = z(inds,i);
        end
    end
    r = rfinal;
    z = zfinal;
else
    if ~rmp.isaxisym
        phi = mod(phi,2*pi/nsym);
        phi(abs(phi - 2*pi/nsym) < 1e-6) = 0;
        inds = find(abs(phi - phistart) < 1e-6);
        r = r(inds,:);
        z = z(inds,:);
    end
end



figure; hold on; box on;
if connect_dots
%     for i = 1:num_surfs
%         plot(r(:,i),z(:,i),'.-')
%     end
    plot(r,z,'.-')
else
    plot(r,z,'.')
end

plot(g.lim(1,:),g.lim(2,:),'k')
plot([g.r(1),g.r(end),g.r(end),g.r(1),g.r(1)],[g.z(1),g.z(1),g.z(end),g.z(end),g.z(1)],'k')
toc