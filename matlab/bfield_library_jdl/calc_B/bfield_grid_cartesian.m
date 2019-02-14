function [Bx,By,Bz]=bfield_grid_cartesian(X,Y,Z,Bgrid,nowarn)
% Bgrid has, X,Y,Z,nx,ny,nz,dX,dY,dZ,Bx,By,Bz, scale_factor
persistent already_warned
if isempty(already_warned)
    already_warned = 0;
end

npts = length(X);
Bx=zeros(npts,1);
By=zeros(npts,1);
Bz=zeros(npts,1);


dX_grid = Bgrid.dX;
dY_grid = Bgrid.dY;
dZ_grid = Bgrid.dZ;

for i = 1:npts
    
    Bx(i) = 0;
    By(i) = 0;
    Bz(i) = 0;   
        
    ix = floor( (X(i) - Bgrid.X(1))/dX_grid ) + 1;
    iy = floor( (Y(i) - Bgrid.Y(1))/dY_grid ) + 1;
    iz = floor( (Z(i) - Bgrid.Z(1))/dZ_grid ) + 1;
    
    if ix < 1 || ix >= Bgrid.nx || iy < 1 || iy >= Bgrid.ny || iz < 1 || iz >= Bgrid.nz
        if ~nowarn && ~already_warned
            warning('Point(s) off grid in bfield_grid_cartesian --> B(:) = 0 there')
            already_warned = 1;
        end
        % already set to zero above
        continue;
    end
    
    z_fac = (Z(i) - Bgrid.Z(iz))/dZ_grid;
    
    dx2 = Bgrid.X(ix+1) - X(i);
    dx1 = dX_grid - dx2;
    dy2 = Bgrid.Y(iy+1) - Y(i);
    dy1 = dY_grid - dy2;    
    
    % Trilinear for each
    QQ1 = Bgrid.Bx(ix:ix+1,iy:iy+1,iz  );
    QQ2 = Bgrid.Bx(ix:ix+1,iy:iy+1,iz+1);
    Bx(i) = ...
        Bgrid.scale_factor*(1-z_fac)*(QQ1(1,1)*dx2*dy2 + QQ1(2,1)*dx1*dy2 + QQ1(1,2)*dx2*dy1 + QQ1(2,2)*dx1*dy1)/(dX_grid*dY_grid) + ...
           Bgrid.scale_factor*z_fac *(QQ2(1,1)*dx2*dy2 + QQ2(2,1)*dx1*dy2 + QQ2(1,2)*dx2*dy1 + QQ2(2,2)*dx1*dy1)/(dX_grid*dY_grid);
    QQ1 = Bgrid.By(ix:ix+1,iy:iy+1,iz  );
    QQ2 = Bgrid.By(ix:ix+1,iy:iy+1,iz+1);
    By(i) = ...
        Bgrid.scale_factor*(1-z_fac)*(QQ1(1,1)*dx2*dy2 + QQ1(2,1)*dx1*dy2 + QQ1(1,2)*dx2*dy1 + QQ1(2,2)*dx1*dy1)/(dX_grid*dY_grid) + ...
           Bgrid.scale_factor*z_fac *(QQ2(1,1)*dx2*dy2 + QQ2(2,1)*dx1*dy2 + QQ2(1,2)*dx2*dy1 + QQ2(2,2)*dx1*dy1)/(dX_grid*dY_grid);
    QQ1 = Bgrid.Bz(ix:ix+1,iy:iy+1,iz  );
    QQ2 = Bgrid.Bz(ix:ix+1,iy:iy+1,iz+1);
    Bz(i) = ...
        Bgrid.scale_factor*(1-z_fac)*(QQ1(1,1)*dx2*dy2 + QQ1(2,1)*dx1*dy2 + QQ1(1,2)*dx2*dy1 + QQ1(2,2)*dx1*dy1)/(dX_grid*dY_grid) + ...
           Bgrid.scale_factor*z_fac *(QQ2(1,1)*dx2*dy2 + QQ2(2,1)*dx1*dy2 + QQ2(1,2)*dx2*dy1 + QQ2(2,2)*dx1*dy1)/(dX_grid*dY_grid);
    
end % npts


