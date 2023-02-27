function P = calc_Bangle_g(g,RZ1,RZ2,npts)
%
% RZ1, RZ2 ... start and end points used to calculate surface normal.
% Surface normal is (RZ2-RZ1,0)x(0,0,1) for (R,Z,phi), with phi coming out
% of the page. This means if wall elements are anticlockwise then the
% surface normal will point "out", away from the axis.
%
% npts is number of points to evaluate along this segment. 
% If npts == 1 the midpoint is used.
%
% P.alpha is the angle of incidence relative to the surface, i.e., the 
% incident flux is qdep = qprl*sin(P.alpha)
%
% P.beta is the angle in the poloidal plane between B and the surface
% normal
%
if length(RZ1) ~= 2 || length(RZ2) ~= 2
    error('Only single point pairs are allowed')
end
if nargin < 3
    npts = 3;
end
if npts < 1
    error('npts must be at least 1')
end

R1 = RZ1(1);
Z1 = RZ1(2);
R2 = RZ2(1);
Z2 = RZ2(2);

DEBUG = 0;
NOWARN = 0;
if DEBUG
    dPlot = 0.01;
    figure; hold on; box on; set(gcf,'color','w');set(gca,'fontsize',14,'fontweight','bold'); grid on;
    axis equal;
    plot([R1,R2],[Z1,Z2],'ko-')
end
v2 = [0,0,1];


if npts == 1
    Reval = 0.5*(R2+R1);
    Zeval = 0.5*(Z2+Z1);
else
    Reval = linspace(R1,R2,npts);
    Zeval = linspace(Z1,Z2,npts);
end

if DEBUG 
    plot(Reval,Zeval,'rx')
end

% Surface normal is the same across the segment
v1 = [Reval(end)-R1,Zeval(end)-Z1,0];
vn = cross(v1,v2);
vn = vn./norm(vn);

% Then local B vector is used
alpha_deg = zeros(1,npts);
beta_deg = zeros(1,npts);
for j = 1:length(Reval)
    if DEBUG
        t = atan2(vn(2),vn(1));
        plot([Reval(j),Reval(j)+dPlot*cos(t)],[Zeval(j),Zeval(j)+dPlot*sin(t)],'b')        
    end
    b = bfield_geq_bicub(g,Reval(j),Zeval(j),NOWARN);
    if isempty(b)
        alpha_deg(j) = NaN;
        beta_deg(j) = NaN;
        continue;
    end
    b = [b.br,b.bz,b.bphi];
       
    b = b./norm(b);
    if DEBUG
        t = atan2(b(2),b(1));
        plot([Reval(j),Reval(j)+dPlot*cos(t)],[Zeval(j),Zeval(j)+dPlot*sin(t)],'m')        
    end


    alpha_deg(j) = asin(dot(vn,b))*180/pi;

    brz = b;
    brz(3) = 0;
    brz = brz./norm(brz);
    beta_deg(j) = acos(dot(vn,brz))*180/pi;
end


P.alpha_deg = alpha_deg;
P.beta_deg = beta_deg;
