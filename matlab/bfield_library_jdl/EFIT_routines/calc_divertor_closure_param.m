clearvars;

gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Qprl experiment\power13mw\g174308.03500_159';
% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\NegD\Proposed\gFile_LSNnegD_fixv2';
% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\NegD\180520-HMode\180520_2800_efit08.gfile';
% gfile_name = 'C:\Work\g185857.04000';
% gfile_name = 'C:\Work\g185857.04000';
% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\EMC3\Florian\g147170.02300';
% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\TDCP\g183299.03000';
% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\TDCP\g184549.02300';
% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\TDCP\g184533.02000';
% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\MAST\mastu_45456_445ms.geqdsk';

iSP_USE = 1

addSpPerp = 0; % Add unit vector perpendicular to Bpol to plots

g = readg_g3d(gfile_name);
psiN_x2 = refine_psi(g.r,g.z,4,g);    % only used for plots
xpt_info = find_xpt_jl(g,1,0,1e-6,1); % only used for plots

%% Make lim (and refined version)
lim = make_lim(g);

%% Find strike points
sp = get_sp(g,lim);

%% Calc some geometric info at the SPs
sp = calc_sp_geo(sp,g,addSpPerp);

%% Plot geometric info
figure; hold on; box on; grid on; set(gcf,'color','w');set(gca,'fontsize',14);
axis equal;
plot(g.lim(1,:),g.lim(2,:),'k','linew',2)
contour(psiN_x2.r,psiN_x2.z,psiN_x2.pn',[1,1],'b-','linewidth',2);
Lplot = .1; % length for drawn unit vectors
for iSP = 1:sp.n
    plot(sp.R(iSP),sp.Z(iSP),'rx','linew',2)
    plot(sp.RSOLSide(iSP),sp.ZSOLSide(iSP),'gx','linew',2)
    plot([sp.R(iSP),Lplot*cos(sp.thetaSurfNorm(iSP))+sp.R(iSP)],[sp.Z(iSP),Lplot*sin(sp.thetaSurfNorm(iSP))+sp.Z(iSP)],'g','linew',2)
    plot([sp.R(iSP),Lplot*cos(sp.thetab(iSP))+sp.R(iSP)],[sp.Z(iSP),Lplot*sin(sp.thetab(iSP))+sp.Z(iSP)],'m','linew',2)
    if addSpPerp
        plot([sp.R(iSP),Lplot*cos(sp.thetabPerpSP(iSP))+sp.R(iSP)],[sp.Z(iSP),Lplot*sin(sp.thetabPerpSP(iSP))+sp.Z(iSP)],'r','linew',2)
        plot([sp.R(iSP),-Lplot*cos(sp.thetabPerpSP(iSP))+sp.R(iSP)],[sp.Z(iSP),-Lplot*sin(sp.thetabPerpSP(iSP))+sp.Z(iSP)],'r','linew',2)
    end
    plot([sp.R(iSP),Lplot*cos(sp.thetaSOL(iSP))+sp.R(iSP)],[sp.Z(iSP),Lplot*sin(sp.thetaSOL(iSP))+sp.Z(iSP)],'y','linew',2)
end
Ltmp = 0.2; % length to pad axis
axis([min([sp.R;xpt_info.rx])-Ltmp,max([sp.R;xpt_info.rx])+Ltmp,min([sp.Z;xpt_info.zx])-Ltmp,max([sp.Z;xpt_info.zx])+Ltmp])


%%
% figure; hold on; box on; grid on; set(gcf,'color','w');set(gca,'fontsize',14);
% contour(psiN_x2.r,psiN_x2.z,psiN_x2.pn',[1,1],'b-','linewidth',2);
% axis equal
% plot(g.lim(1,:),g.lim(2,:),'k-','linewidth',2)
% plot(sp.R,sp.Z,'rx','linew',2)




%%
LminT = 0.5;
nRays = 100;
theta_testOSP = linspace(1e-3,sp.alpha_deg_polplane(iSP_USE)*pi/180,nRays);
[LtoIntOSP,PintOSP] = check_Rays(theta_testOSP, sp, iSP_USE, g, lim, 1);

T = ones(size(LtoIntOSP));
T(LtoIntOSP > LminT) = 0;
% FFun = @(theta) sin(theta).^2;
FFun = @(theta) sin(theta);
% FFun = @(theta) ones(size(theta));
F = FFun(theta_testOSP);
tInt = linspace(0,pi,1000);
FInt = FFun(tInt);
FNorm = trapz(tInt,FInt);
plot(PintOSP(T == 1,1),PintOSP(T == 1,2),'go')


% figure; hold on; box on; grid on; set(gcf,'color','w');set(gca,'fontsize',14);
% axis equal;
% plot(g.lim(1,:),g.lim(2,:),'k','linew',2)
% contour(psiN_x2.r,psiN_x2.z,psiN_x2.pn',[1,1],'b-','linewidth',2);

% theta_testFull = linspace(1e-3,pi-1e-3,nRays);
% [LtoIntFull,PintFull] = check_Rays(theta_testFull, sp, iSP_USE, g, lim, 1);


% TFull = ones(size(LtoIntFull));
% TFull(LtoIntFull > LminT) = 0;
% FFull = sin(theta_testFull).^2;
% % FFull = sin(theta_testFull*2);
% plot(PintFull(TFull == 1,1),PintFull(TFull == 1,2),'go')

%%
Tint = trapz(theta_testOSP,T)/max(theta_testOSP); % geometric fraction of particles going towards div < Lmin
Fint = trapz(theta_testOSP,F)/FNorm; % Fraction of cosine distribution particles on SOL side
Totint = trapz(theta_testOSP,F.*T)/trapz(theta_testOSP,F);

% TintFull = trapz(theta_testFull,TFull)/max(theta_testFull); % geometric fraction of particles going towards div < Lmin
% FintFull = trapz(theta_testFull,FFull)/FNorm; % Fraction of cosine distribution particles on SOL side
% TotintFull = trapz(theta_testFull,FFull.*TFull)/trapz(theta_testFull,FFull);


fprintf('T integrated (Lmin = %.1f) : %.5f\n',LminT,Tint)
fprintf('F integrated              : %.5f\n',Fint)
fprintf('F*T integrated            : %.5f\n',Totint)

% fprintf('full:\n')
% fprintf('T integrated (Lmin = %.1f) : %.5f\n',LminT,TintFull)
% fprintf('F integrated              : %.5f\n',FintFull)
% fprintf('F*T integrated            : %.5f\n',TotintFull)


% fprintf('cos(beta)*F*T int         : %.5f\n',sind(sp.beta_deg(iSP_USE))*Totint)

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function sp = calc_sp_geo(sp,g,addSpPerp)
for iSP = 1:sp.n
    % Provide a point on the SOL side: check both lim pts and find the one
    % with the same sign
    ptest = calc_psiN(g,[sp.R(iSP),sp.R1Fine(iSP),sp.R2Fine(iSP)],[sp.Z(iSP),sp.Z1Fine(iSP),sp.Z2Fine(iSP)]);
    sp.iUse(iSP) = find(ptest(2:end) > 1);

    if sp.iUse(iSP) == 1
        sp.RSOLSide(iSP) = sp.R1Fine(iSP);
        sp.ZSOLSide(iSP) = sp.Z1Fine(iSP);
    else
        sp.RSOLSide(iSP) = sp.R2Fine(iSP);
        sp.ZSOLSide(iSP) = sp.Z2Fine(iSP);
    end

    % make unit vectore along PFC towards SOL side
    v = [sp.RSOLSide(iSP)-sp.R(iSP),sp.ZSOLSide(iSP)-sp.Z(iSP),0];
    v = v./norm(v);
    sp.vSOL(iSP,:) = v;
    sp.thetaSOL(iSP) = atan2(sp.vSOL(iSP,2),sp.vSOL(iSP,1));


    % Get surface normal and bfield normal
    v2 = [0,0,1];
    v1 = [sp.R(iSP)-sp.RSOLSide(iSP),sp.Z(iSP)-sp.ZSOLSide(iSP),0];
    vn = cross(v1,v2);
    vn = vn./norm(vn);

    b = bfield_geq_bicub(g,sp.R(iSP),sp.Z(iSP),0);
    b = [b.br,b.bz,b.bphi];
    b = b./norm(b);

    if addSpPerp
        % normal to b at SP
        bPerpSP = cross(b,[0,0,1]);
        sp.bPerpSP(iSP,:) = bPerpSP./norm(bPerpSP);
    end
    sp.alpha_deg(iSP) = 90 - acos(dot(vn,b))*180/pi; % Actual qperp = q||*sin(alpha)
    sp.vSurfNorm(iSP,:) = vn;
    sp.bTot(iSP,:) = b;

    %% local angles of unit vectors
    sp.thetaSurfNorm(iSP) = atan2(vn(2),vn(1));
    sp.thetab(iSP) = atan2(b(2),b(1));
    if addSpPerp
        sp.thetabPerpSP(iSP) = atan2(sp.bPerpSP(iSP,2),sp.bPerpSP(iSP,1));
    end
    %% check sign of surface norm
    % Take a small step along norm and check if inside vessel polygon
    Lstep = 1e-3;
    pCheck = [Lstep.*cos(sp.thetaSurfNorm(iSP))+sp.R(iSP),Lstep.*sin(sp.thetaSurfNorm(iSP))+sp.Z(iSP)];
    isIn = inpolygon(pCheck(1),pCheck(2),g.lim(1,:),g.lim(2,:));
    if ~isIn
        sp.vSurfNorm(iSP,:) = -sp.vSurfNorm(iSP,:);
        sp.thetaSurfNorm(iSP) = atan2(sp.vSurfNorm(iSP,2),sp.vSurfNorm(iSP,1));
    end

    brz = b;
    brz(3) = 0;
    sp.bPol(iSP,:) = brz./norm(brz);

    %% Change sign on bPol to point towards x-point if necessary
    % Take a small step along norm and check if inside vessel polygon
    Lstep = 1e-3;
    pCheck = [Lstep.*cos(sp.thetab(iSP))+sp.R(iSP),Lstep.*sin(sp.thetab(iSP))+sp.Z(iSP)];
    isIn = inpolygon(pCheck(1),pCheck(2),g.lim(1,:),g.lim(2,:));
    if ~isIn
        sp.bPol(iSP,:) = -sp.bPol(iSP,:);
        sp.thetab(iSP) = atan2(sp.bPol(iSP,2),sp.bPol(iSP,1));
    end

    if addSpPerp
        sp.beta_deg(iSP) = acosd(dot(sp.vSOL(iSP,:),sp.bPerpSP(iSP,:)));
        if sp.beta_deg(iSP) > 90
            sp.beta_deg(iSP) = 180 - sp.beta_deg(iSP);
        end
    end
    sp.alpha_deg_polplane(iSP) = acosd(dot(sp.vSOL(iSP,:),sp.bPol(iSP,:)));
    % sp.pit.*sind(sp.beta_deg)  should be sind(sp.alpha_deg)

end

sp.B = bfield_geq_bicub(g,sp.R,sp.Z);
sp.B.bpol = sqrt(sp.B.bz.^2 + sp.B.br.^2);
sp.pit = atan2(sp.B.bpol,abs(sp.B.bphi));
sp.pit_deg = atan2d(sp.B.bpol,abs(sp.B.bphi));

end


%%
function sp = get_sp(g,lim)
[pint,~,~,~,int_count] = int_curve_curve(lim.Lfine,lim.psiNfine,[0,max(lim.L)],[1,1]);
if ~any(int_count == [2,4])
    plot_gfile(g)
    Lint = pint(:,1);
    [~,R,Z,iLim] = psiN_at_L(g,lim.L,lim.R,lim.Z,Lint);
    plot(R(1),Z(1),'gx','linew',4)
    plot(R(2:end),Z(2:end),'co')
    error('Why not two or four intersections?')
end
sp.n = int_count;
Lint = pint(:,1);

% These are the indices of the g lim structure on either side of the SP
[~,sp.R,sp.Z,iLim] = psiN_at_L(g,lim.L,lim.R,lim.Z,Lint);
iLim1 = floor(iLim);
iLim2 = ceil(iLim);
sp.R1 = lim.R(iLim1);
sp.Z1 = lim.Z(iLim1);
sp.R2 = lim.R(iLim2);
sp.Z2 = lim.Z(iLim2);

[~,~,~,iLimFine] = psiN_at_L(g,lim.Lfine,lim.Rfine,lim.Zfine,Lint);
iLim1Fine = floor(iLimFine);
iLim2Fine = ceil(iLimFine);
sp.R1Fine = lim.Rfine(iLim1Fine);
sp.Z1Fine = lim.Zfine(iLim1Fine);
sp.R2Fine = lim.Rfine(iLim2Fine);
sp.Z2Fine = lim.Zfine(iLim2Fine);
end




%%
% figure; hold on;
% % plot(lim.L,lim.psiN,'o-')
% plot(lim.Lfine,lim.psiNfine,'x-')
% yline(1)
% xline(pint1(:,1))

%%
function lim = make_lim(g)
lim.R = g.lim(1,1:end);
lim.Z = g.lim(2,1:end);
lim.dL = sqrt(diff(lim.R).^2 + diff(lim.Z).^2);
tol = 1e-6;
iTol = find(lim.dL < tol);
lim.R(iTol) = []; lim.Z(iTol) = []; lim.dL(iTol) = [];
lim.L = [0,cumsum(lim.dL)];
lim.psiN = calc_psiN(g,lim.R,lim.Z)';

lim.Lfine = linspace(0,max(lim.L),1000);
% Need to make sure actual lim vertex points are in there, otherwise the
% corners are smoothed!
lim.Lfine = sort([lim.Lfine,lim.L(2:end-1)]);
iTol = find(diff(lim.Lfine) < tol);
lim.L(iTol) = [];
[lim.psiNfine,lim.Rfine,lim.Zfine] = psiN_at_L(g,lim.L,lim.R,lim.Z,lim.Lfine);
end


%%
function [psiN,r,z,i] = psiN_at_L(g,limL,limR,limZ,Lwant)
r = interp1(limL,limR,Lwant);
z = interp1(limL,limZ,Lwant);
i = interp1(limL,1:length(limL),Lwant);
psiN = calc_psiN(g,r,z);
end

%
function [psiN_x2] = refine_psi(r,z,fac,g)
psiN_x2.r = linspace(r(2),r(end-1),length(r)*fac);
psiN_x2.z = linspace(z(2),z(end-1),length(z)*fac);
for i = 1:length(psiN_x2.z)
    psiN_x2.pn(:,i) = calc_psiN(g,psiN_x2.r,psiN_x2.z(i)*ones(size(psiN_x2.r)));
end
end

function [LtoInt,pInts] = check_Rays(theta_test, sp, iSP_USE, g, lim,plotRays)


% plotRays = 1;
Ltest = 3;

rSP = sp.R(iSP_USE);
zSP = sp.Z(iSP_USE);


% Calculate the angle from the x-axis to the line tang to the surface
% so need a point in the right direction along the PFC
if sp.iUse(iSP_USE) == 1
    rPFC = sp.R1(iSP_USE);
    zPFC = sp.Z1(iSP_USE);
else
    rPFC = sp.R2(iSP_USE);
    zPFC = sp.Z2(iSP_USE);
end
gamma = atan2(zPFC-zSP,rPFC-rSP);
gamma = mod(gamma + 2*pi,2*pi);

for i = 1:length(theta_test)
    % Figure out what direction to go
    angTotCCW = gamma + theta_test(i);
    angTotCW = gamma - theta_test(i);
    LCWtest = 1e-3;
    rCWtest = LCWtest*cos(angTotCW) + rSP;
    zCWtest = LCWtest*sin(angTotCW) + zSP;
    rCCWtest = LCWtest*cos(angTotCCW) + rSP;
    zCCWtest = LCWtest*sin(angTotCCW) + zSP;
    isInCCW = inpolygon(rCCWtest,zCCWtest,g.lim(1,:),g.lim(2,:));
    isInCW = inpolygon(rCWtest,zCWtest,g.lim(1,:),g.lim(2,:));
    if ~xor(isInCCW,isInCW)
        error('Must be exclusive')
    end
    if isInCCW
        angTot = angTotCCW;
    else
        angTot = angTotCW;
    end
    % end point of ray
    r4 = Ltest*cos(angTot) + rSP;
    z4 = Ltest*sin(angTot) + zSP;

    %     if i == 2
    %         plot([rSP,r4],[zSP,z4],'-','linew',2)
    %     end

    % Move a small distance along the test line to avoid self-intersection
    rSPoff = 0.001*cos(angTot) + rSP;
    zSPoff = 0.001*sin(angTot) + zSP;

    % Find FIRST intersection with PFC contour
    [pint1,ierr,found_ind,int_count]=int_line_curve([rSPoff,zSPoff],[r4,z4],lim.R,lim.Z);
    if isnan(pint1(1))
        plot([rSP,r4],[zSP,z4],'k-','linew',2)
        error('increase Ltest')
    end
    if int_count ~=1
        LtoInt_tmp = sqrt((rSP-pint1(:,1)).^2 + (zSP-pint1(:,2)).^2);
        [LtoInt(i),imin] = min(LtoInt_tmp);
        pint1(1,:) = pint1(imin,:);
    else
        LtoInt(i) = sqrt((rSP-pint1(1,1))^2 + (zSP-pint1(1,2))^2);
    end
    pInts(i,:) = pint1(1,1:2);
    if plotRays
        plot(pint1(1,1),pint1(1,2),'bo','linew',1)
        plot([rSP,pint1(1,1)],[zSP,pint1(1,2)],'k-')
    end
end
end