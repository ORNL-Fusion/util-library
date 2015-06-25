function [izc,irc,jpc,ierr]=point2cell_c2(grid,r0,z0,quiet)
% RZ in meters

if nargin < 4
    quiet = 0;
end

debug_text = 0;
debug_plot = 0;
npts = length(r0);
ierr = zeros(1,npts);

for ipt = 1:npts
    
    % Identify zone by checking in polygon with boundary
    inpoly0 = zeros(1,grid.ndomain);
%     figure; hold on;
    for iz = 1:grid.ndomain
        nr = grid.ny(iz);
        np = grid.nx(iz);
        
        rpoly0 = [grid.x2d{iz}(1:nr,1);grid.x2d{iz}(nr,2:np).';grid.x2d{iz}(nr-1:-1:1,np);grid.x2d{iz}(1,np-1:-1:1).'];
        zpoly0 = [grid.y2d{iz}(1:nr,1);grid.y2d{iz}(nr,2:np).';grid.y2d{iz}(nr-1:-1:1,np);grid.y2d{iz}(1,np-1:-1:1).'];        
        inpoly0(iz) = inpolygon(r0(ipt),z0(ipt),rpoly0,zpoly0);
        
%             plot(rpoly0,zpoly0)
%             plot(r0(ipt),z0(ipt),'rx')
    end
    
    if sum(inpoly0) ~= 1
        if sum(inpoly0) == 0 && ~quiet
            fprintf('Point not found in any zone!\n')
        end
        if sum(inpoly0) > 1
            disp(['Point found in more than one zone!? : ',num2str(inpoly0)])
        end
        
        ierr(ipt) = 1;
        izc(ipt)=NaN;irc(ipt)=NaN;jpc(ipt)=NaN;
        continue;
        %     return
    end
    
    %%%%%%%%%%%
    izc(ipt) = find(inpoly0 == 1);
    %%%%%%%%%%%
    
    % Identify radial and poloidal cell indices
    nr = grid.ny(izc(ipt));
    np = grid.nx(izc(ipt));
    foundit = 0;
    for ir = 1:nr-1
        for jp = 1:np-1
            rpoly0 = [grid.x2d{izc(ipt)}(ir,jp,ipt),grid.x2d{izc(ipt)}(ir+1,jp,ipt),grid.x2d{izc(ipt)}(ir+1,jp+1,ipt),grid.x2d{izc(ipt)}(ir,jp+1,ipt),grid.x2d{izc(ipt)}(ir,jp,ipt)];
            zpoly0 = [grid.y2d{izc(ipt)}(ir,jp,ipt),grid.y2d{izc(ipt)}(ir+1,jp,ipt),grid.y2d{izc(ipt)}(ir+1,jp+1,ipt),grid.y2d{izc(ipt)}(ir,jp+1,ipt),grid.y2d{izc(ipt)}(ir,jp,ipt)];
            inpoly0 = inpolygon(r0(ipt),z0(ipt),rpoly0,zpoly0);
            
            if inpoly0 == 1
                foundit = 1;
                irc(ipt) = ir;
                jpc(ipt) = jp;
                break;
            end
        end
        if foundit == 1
            break
        end
    end
    
    if foundit == 0
        ierr(ipt) = 1;
        fprintf('Point not found in any cell!\n')
        izc(ipt)=NaN;irc(ipt)=NaN;jpc(ipt)=NaN;
        continue;
    end
    
    % Debug 
    if debug_text        
        fprintf('Found point [R,Z] = [%8.3f, %8.3f] at [iz,ir,jp] = [%d,%d,%d]\n',[r0(ipt),z0(ipt),izc(ipt),irc(ipt),jpc(ipt)])
    end
    if debug_plot
        ZNP = izc(ipt);
        IR = irc(ipt);
        JP = jpc(ipt);        
        
        rpoly0 = [grid.x2d{ZNP}(IR,JP),grid.x2d{ZNP}(IR+1,JP),grid.x2d{ZNP}(IR+1,JP+1),grid.x2d{ZNP}(IR,JP+1),grid.x2d{ZNP}(IR,JP)];
        zpoly0 = [grid.y2d{ZNP}(IR,JP),grid.y2d{ZNP}(IR+1,JP),grid.y2d{ZNP}(IR+1,JP+1),grid.y2d{ZNP}(IR,JP+1),grid.y2d{ZNP}(IR,JP)];
        figure; hold on;        
        plot(rpoly0,zpoly0,'.k--')
        plot(r0(ipt),z0(ipt),'rx')
    end        
end