function [coil,COIL,coil_emc3] = make_larry_coils(nfil_coil,R_offset)
%
%  Note: lower case r,z,phi,x,y is global cyclindrical system, 'machine coordinates', 
%  where z is vertical through a tokamak axis of symmetry;
%  Upper case R,Z,PHI,X,Y has Z running down the linear machine axis

if nargin < 2
    R_offset = 0;
end

% nfil_coil = 10;

plot_schematic = 0;
plot_coils = 0;

% LARRY DATA -  LARRY DATA -  LARRY DATA -  LARRY DATA -  LARRY DATA -  LARRY DATA
NCOIL = 15;
R1 = [  0.32,  0.075,   0.32,  0.075,   0.32,  0.075,   0.32,  0.075,   0.32,  0.075,   0.32,  0.075,   0.32,  0.075,   0.32];
Z1 = [-1.005, -0.505, -0.005,  0.495,  0.995,  1.495,  1.995,  2.495,  2.995,  3.495,  3.995,  4.495,  4.995,  5.495,  5.995];
R2 = [  0.33,  0.085,   0.33,  0.085,   0.33,  0.085,   0.33,  0.085,   0.33,  0.085,   0.33,  0.085,   0.33,  0.085,   0.33];
Z2 = [-0.995, -0.495,  0.005,  0.505,  1.005,  1.505,  2.005,  2.505,  3.005,  3.505,  4.005,  4.505,  5.005,  5.505,  6.005];
CURRENT = 4000;
WINDINGS = 4;
% LARRY DATA -  LARRY DATA -  LARRY DATA -  LARRY DATA -  LARRY DATA -  LARRY DATA

curr = CURRENT*WINDINGS;

%
% Plot Larry's coils, and define radius and Z for circular filaments
%
if plot_schematic == 1
    figure; hold on; box on;
end
for ic = 1:NCOIL
    cR = [R1(ic), R2(ic), R2(ic), R1(ic), R1(ic)];
    cZ = [Z1(ic), Z1(ic), Z2(ic), Z2(ic), Z1(ic)];
    if plot_schematic == 1
        plot(cR,cZ,'k','linewidth',2)
    end
    cR_flip = -[R1(ic), R2(ic), R2(ic), R1(ic), R1(ic)];
    cZ_flip = [Z1(ic), Z1(ic), Z2(ic), Z2(ic), Z1(ic)];    
    
    if plot_schematic == 1
        plot(cR_flip,cZ_flip,'b','linewidth',2)
    end
    Rcen(ic) = mean(cR(1:4));
    Zcen(ic) = mean(cZ(1:4));
    if plot_schematic == 1
        plot(Rcen(ic),Zcen(ic),'r.')
    end
end
Rmax = max([R1,R2]);
Zmax = max([Z1,Z2]);
Zmin = min([Z1,Z2]);
if plot_schematic == 1
    axis([-Rmax*1.1,Rmax*1.1,Zmin*1.1,Zmax*1.1])
    plot([0.07,0.07,-0.07,-0.07,0.07],[0,5,5,0,0],'k--','linewidth',2)
    xlabel('R (m)','fontsize',12)
    ylabel('Z (m)','fontsize',12)
end
%
% Make coils out of filaments
%
if plot_coils == 1
figure; hold on; box on;
end
for ic = 1:NCOIL
    Rcoil = ones(1,nfil_coil+1)*Rcen(ic);
    Zcoil = ones(1,nfil_coil+1)*Zcen(ic);
    Pcoil = linspace(0,2*pi,nfil_coil+1);    
    Xcoil = Rcoil.*cos(Pcoil);
    Ycoil = Rcoil.*sin(Pcoil);
           
    % convert to machine coords
    xcoil = Rcoil.*cos(Pcoil);
    zcoil = Rcoil.*sin(Pcoil);
    ycoil = Zcoil;
    
%     rcoil = sqrt(xcoil.^2 + ycoil.^2);
%     pcoil = atan2(ycoil,xcoil);
    
    if plot_coils == 1
    plot3(xcoil,ycoil,zcoil,'r.-')
    end
    inds = 1+(ic-1)*(nfil_coil+1):(ic)*(nfil_coil+1);
    coil.coil(inds,1) = xcoil;
    coil.coil(inds,2) = ycoil;
    coil.coil(inds,3) = zcoil;
    coil.current(inds(1:end-1),1) = curr;
    coil.current(inds(end)) = 0;
    
    COIL.coil(inds,1) = Xcoil;
    COIL.coil(inds,2) = Ycoil;
    COIL.coil(inds,3) = Zcoil;    
    
%     c.rcoil(inds) = rcoil;
%     c.zcoil(inds) = zcoil;
end
COIL.current = coil.current;
if plot_coils == 1
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
end


coil_emc3.current = coil.current;

xx = COIL.coil(:,1);
yy = COIL.coil(:,2);
zz = COIL.coil(:,3)-2.5;
% xx = coil.coil(:,1);
% yy = coil.coil(:,2)-2.5;
% zz = coil.coil(:,3);
th = atan2(yy,xx);
rr = sqrt(xx.^2 + yy.^2);

rrr = R_offset + rr.*cos(th);
zzz = rr.*sin(th);
ppp = atan2(zz,rrr);
xxx = rrr.*cos(ppp);
yyy = rrr.*sin(ppp);

coil_emc3.coil(:,1) = xxx;
coil_emc3.coil(:,2) = yyy;
coil_emc3.coil(:,3) = zzz;

% figure; hold on; box on;
% plot3(coil.coil(:,1),coil.coil(:,2),coil.coil(:,3),'k.-')
% xlabel('x (m)')
% ylabel('y (m)')
% zlabel('z (m)')
% 
% figure; hold on;
% plot(c.rcoil,c.zcoil)


% plot3(COIL.coil(:,1),COIL.coil(:,2),COIL.coil(:,3),'k.-')
    