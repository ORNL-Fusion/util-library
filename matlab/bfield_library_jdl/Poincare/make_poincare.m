function poinc = make_poincare(bfield,Rstart,Zstart,phistart,phi_want,npoints_want,dphi,sort_it,Rax_guess,Zax_guess,plot_settings,save_name)
% R,Z = [m], phi = [radians]
% R,Z can be arrays
if nargin < 5
    phi_want = phistart;
end
if nargin < 6
    npoints_want = 10;
end
if nargin < 7
    dphi = 0.5*pi/180;
end
if nargin < 8
    sort_it = 0;
end
if nargin >= 12
    saveit = 1;
else
    saveit = 0;
end

if saveit
    if exist(save_name,'file') == 2
        inp = input('The output file already exists, are you sure you want to overwrite [Y/N] ??','s');
        if ~strcmpi(inp,'y')
            error('Stopping')
        end
    end
end


tic;
fprintf('Making %d poincare surfaces\n',length(Rstart))

%--------------------------------------------------------------------------
%      POINCARE SETUP
%--------------------------------------------------------------------------
nsurfs = length(Rstart);
if nsurfs ~= length(Zstart)
    if length(Zstart) == 1
        Zstart = Zstart*ones(size(Rstart));
    else
        error('Input inconsistency')
    end
end
if length(phistart) ~= 1
    error('No')
end

NFP = bfield.nsym;
ntran = npoints_want/NFP;
nsteps_real = abs(ntran*2*pi/dphi);
nsteps = abs(round(nsteps_real));
if abs(nsteps_real - nsteps) > 1e-6
    error('Make this divisible')
end

nowarn = 0;
istart = (phi_want - phistart)/dphi + 1;
istride = 2*pi/NFP/dphi;
npoinc = length(istart:istride:nsteps+1);
Rpoinc = NaN(npoinc,nsurfs);
Zpoinc = NaN(npoinc,nsurfs);
for i = 1:nsurfs
    fprintf('Working on surface %d of %d\n',i,nsurfs);
    f = follow_fieldlines_rzphi_dphi(bfield,Rstart(i),Zstart(i),phistart,dphi,nsteps,nowarn);
    Rpoinc(:,i) = f.r(istart:istride:end,1);
    Zpoinc(:,i) = f.z(istart:istride:end,1);   
    if 0
        x = f.r.*cos(f.phi);
        y = f.r.*sin(f.phi);
        z = f.z;
        figure; hold on; box on;
        plot3(x,y,z)
    end
end

fprintf('Poincare took %f seconds.\n',toc); tic;

%--------------------------------------------------------------------------
%      PLOTTING
%--------------------------------------------------------------------------


if sort_it
%     [Rax,Zax] = find_axis(bfield,phistart,Rax_guess,Zax_guess);
    warning('NOT FINDING AXIS FOR SORTING, USING GUESS!!')
    Rax = Rax_guess;
    Zax = Zax_guess;
    tpoinc = atan2(Zpoinc-Zax,Rpoinc-Rax);
    rpoinc = sqrt((Rpoinc-Rax).^2 + (Zpoinc-Zax).^2);
    [tpoinc,Isort]=sort(tpoinc);                       %sort by increasing theta
    for j = 1:nsurfs
        rpoinc(:,j) = rpoinc(Isort(:,j),j);
        Rpoinc(:,j) = Rpoinc(Isort(:,j),j);
        Zpoinc(:,j) = Zpoinc(Isort(:,j),j);
    end
end
    
if plot_settings.plotit
    if plot_settings.newfig
    figure; hold on; box on;
    end    
    if sort_it
        if plot_settings.connect
            plot(Rpoinc,Zpoinc,'.-','linewidth',2)        
        else
            plot(Rpoinc,Zpoinc,'.','linewidth',2)
        end
    else
        plot(Rpoinc,Zpoinc,'o')
    end
    xlabel('R [m]','fontsize',14)
    ylabel('Z [m]','fontsize',14)
    set(gca,'fontsize',14);
    title(strcat('\phi = ',sprintf('%f',phi_want*180/pi)));
end
poinc.Rpoinc = Rpoinc;
poinc.Zpoinc = Zpoinc;
poinc.phi = phi_want;

save(save_name,'poinc')