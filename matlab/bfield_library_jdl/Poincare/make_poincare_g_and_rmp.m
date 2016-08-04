clearvars;
tic;
connect_dots = 0;

gfile_name = 'C:\Work\DIII-D\164723\g164723.03059_410';
g = readg_g3d(gfile_name);

% bfield.type = 'gfile';
% bfield.g = g;
% isaxisym = 1;

bfield.g = g;
% rmp = build_d3d_icoils_jl([-2903.   2939.  -2889.   2935.  -2886.   2940. -2851.   2907.  -2866.   2918.  -2910.   2918.]);
taper = [-2900.d0 2900.d0 -2900.d0 2900.d0 -2900.d0 2900.d0 -2900.d0 2900.d0 -2900.d0 2900.d0 -2900.d0 2900.d0];
rmp = build_d3d_icoils_jl(taper);

bfield.type = 'gfile+coils';
bfield.coil = rmp.coil;
bfield.current = rmp.current;
isaxisym = 0;
nsym = 3;


Rend   = 2.25;    % if [] this is set to max radius of g.bdry
Rstart = 2.05;    % if [] this is set to g.rmaxis

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


s1=follow_fieldlines_rzphi_dphi(bfield,Rstarts,Zstarts,phistart,dphi,nsteps);
s2=follow_fieldlines_rzphi_dphi(bfield,Rstarts,Zstarts,phistart,-dphi,nsteps);
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
        if ~isaxisym
            phi = mod(phi,2*pi/nsym);
            inds = find(abs(phi - phistart_deg) < 1e-6);
            rfinal(:,i) = r(inds,i);
            zfinal(:,i) = z(inds,i);
        end
    end
    r = rfinal;
    z = zfinal;
else
    if ~isaxisym
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