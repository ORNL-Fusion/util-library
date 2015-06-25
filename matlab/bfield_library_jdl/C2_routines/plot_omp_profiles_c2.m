function plot_omp_profiles_c2(grid,data)
% Currently set up for d3d case and set to give max 1 point per cell
rax = 2; rend = 2.5; zax = 0; zend = 0; 

nr_prof = 200;

data_plot1 = data.Te2d; fac1 = 1.e-3;
data_plot2 = data.Ti2d; fac2 = 1.e-3;
data_plot3 = data.np2d; fac3 = 1;

% intersect line with core inner surface and SOL outer surface
rmin = 100;
rmax = 0;
for iz = 1:grid.ndomain
    [pint1,ierr1]=int_line_curve([rax,zax],[rend,zend],grid.x2d{iz}(1,:),grid.y2d{iz}(1,:));
    if ierr1 == 0
        rmin = min([rmin,pint1(1)]);
        rmax = max([rmax,pint1(1)]);
    end
    [pint1,ierr1]=int_line_curve([rax,zax],[rend,zend],grid.x2d{iz}(grid.ny(iz),:),grid.y2d{iz}(grid.ny(iz),:));
    if ierr1 == 0
        rmin = min([rmin,pint1(1)]);
        rmax = max([rmax,pint1(1)]);
    end    
end

r_prof = linspace(rmin-1e-3,rmax+1e-3,nr_prof);
z_prof = linspace(zax,zax,nr_prof);

data_prof1 = NaN(1,nr_prof);
data_prof2 = NaN(1,nr_prof);
data_prof3 = NaN(1,nr_prof);

ig_used = []; icount = 0;
for irp = 1:nr_prof
    [izc,irc,jpc,ierrc]=point2cell_c2(grid,r_prof(irp),z_prof(irp),1);
    if ierrc == 0
        ig = irc + (jpc-1)*(grid.ny(izc)-1) + grid.Zoffsets(izc);        
        data_prof1(irp) = data_plot1{izc}(irc,jpc);
        data_prof2(irp) = data_plot2{izc}(irc,jpc);
        data_prof3(irp) = data_plot3{izc}(irc,jpc);
        
        if ~any(ig_used == ig)
            icount = icount + 1;
            r_prof_1pc(icount) = r_prof(irp);
            z_prof_1pc(icount) = z_prof(irp);
            data_prof1_1pc(icount) = data_plot1{izc}(irc,jpc);
            data_prof2_1pc(icount) = data_plot2{izc}(irc,jpc);
            data_prof3_1pc(icount) = data_plot3{izc}(irc,jpc);
        end
        ig_used = [ig_used,ig];
        
    end
end


figure; hold on; box on;
plot(r_prof_1pc,data_prof1_1pc*fac1,'ro-','linewidth',2)
plot(r_prof_1pc,data_prof2_1pc*fac2,'bo-','linewidth',2)
plot(r_prof_1pc,data_prof3_1pc*fac3,'ko-','linewidth',2)
xlabel('R_{OMP} [m]','fontsize',12)

ax_use = axis(gca); ax_use(3) = 0; axis(ax_use);

ylabel('T (keV), n (10^{19} m^{-3})','fontsize',12)
legend('T_e','T_i','n_e');
set(gca,'fontsize',12)
