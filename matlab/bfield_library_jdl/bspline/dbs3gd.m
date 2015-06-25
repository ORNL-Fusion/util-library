function val = dbs3gd(iderx,idery,iderz,nxvec,xvec,nyvec,yvec,nzvec,zvec,kx,ky,kz,xknot,yknot,zknot,nx,ny,nz,bcoef,ldf,mdf)
%
% %
% %  Evaluates the derivative of a three-dimensional tensor-product spline,
% %  given its tensor-product B-spline representation on a grid.
% %
% %   iderx  - order of the x-derivative.  (input)
% %   idery  - order of the y-derivative.  (input)
% %   iderz  - order of the z-derivative.  (input)
% %   nx     - number of grid points in the x-direction.  (input)
% %   xvec   - array of length nx containing the x-coordinates at
% %            which the spline is to be evaluated.  (input)
% %            the points in xvec should be strictly increasing.
% %   ny     - number of grid points in the y-direction.  (input)
% %   yvec   - array of length ny containing the y-coordinates at
% %            which the spline is to be evaluated.  (input)
% %            the points in yvec should be strictly increasing.
% %   nz     - number of grid points in the z-direction.  (input)
% %   zvec   - array of length nz containing the z-coordinates at
% %            which the spline is to be evaluated.  (input)
% %            the points in yvec should be strictly increasing.
% %   kx     - order of the spline in the x-direction.  (input)
% %   ky     - order of the spline in the y-direction.  (input)
% %   kz     - order of the spline in the z-direction.  (input)
% %   xknot  - array of length nx+kx containing the knot
% %            sequence in the x-direction.  (input)
% %            xknot must be nondecreasing.
% %   yknot  - array of length ny+ky containing the knot
% %            sequence in the y-direction.  (input)
% %            yknot must be nondecreasing.
% %   zknot  - array of length nz+kz containing the knot
% %            sequence in the z-direction.  (input)
% %            zknot must be nondecreasing.
% %   nx     - number of B-spline coefficients in the x-direction.
% %            (input)
% %   ny     - number of B-spline coefficients in the y-direction.
% %            (input)
% %   nz     - number of B-spline coefficients in the z-direction.
% %            (input)
% %   bcoef  - array of length nx*ny*nz containing the
% %            tensor-product B-spline coefficients.  (input)
% %            bscoef is treated internally as a matrix of size nx
% %            by ny by nz.
% %   val    - array of size nx by ny by nz containing the values of
% %            the (iderx,idery,iderz) derivative of the spline on the
% %            nx by ny by nz grid.  (output)
% %            value(i,j,k) contains the derivative of the spline at
% %            the point (xvec(i), yvec(j), zvec(k)).
% %   ldf    - leading dimension of value exactly as specified in the
% %            dimension statement of the calling program.  (input)
% %   mdf    - middle dimension of value exactly as specified in the
% %            dimension statement of the calling program.  (input)
% %
%
%     integer, intent(in)                               :: iderx, idery, iderz
%     integer, intent(in)                               :: nxvec, nyvec, nzvec
%     integer, intent(in)                               :: kx, nx, ky, ny, kz, nz
%     integer, intent(in)                               :: ldf,mdf
%
%     real(kind=dbl), dimension(nxvec), intent(in)      :: xvec
%     real(kind=dbl), dimension(nyvec), intent(in)      :: yvec
%     real(kind=dbl), dimension(nzvec), intent(in)      :: zvec
%     real(kind=dbl), dimension(nx+kx), intent(in)      :: xknot
%     real(kind=dbl), dimension(ny+ky), intent(in)      :: yknot
%     real(kind=dbl), dimension(nz+kz), intent(in)      :: zknot
%     real(kind=dbl), dimension(nx,ny,nz), intent(in)   :: bcoef
%     real(kind=dbl), dimension(ldf,mdf,*), intent(out) :: val
%
%     integer                                           :: i, ik, il, ix, iy, iz
%     integer                                           :: ikx, iky, ikz
%     integer, dimension(nxvec)                         :: leftx
%     integer, dimension(nyvec)                         :: lefty
%     integer, dimension(nzvec)                         :: leftz
%     real(kind=dbl), dimension(nxvec,kx)               :: biatx
%     real(kind=dbl), dimension(nyvec,ky)               :: biaty
%     real(kind=dbl), dimension(nzvec,kz)               :: biatz
%     real(kind=dbl), dimension(max(nxvec,nyvec,nzvec)) :: term, save1
%
%     real(kind=dbl), dimension(max(nxvec,nyvec,nzvec), max(kx,ky,kz)) :: dl, dr
%
%     logical :: same,next

leftx = zeros(1,nxvec);
lefty = zeros(1,nyvec);
leftz = zeros(1,nzvec);
biatx = zeros(nxvec,kx);
biaty = zeros(nyvec,ky);
biatz = zeros(nzvec,kz);
save1 = zeros(1,max([nxvec,nyvec,nzvec]));
term = zeros(1,max([nxvec,nyvec,nzvec]));
dl = zeros(max([nxvec,nyvec,nzvec]), max([kx,ky,kz]));
dr = zeros(max([nxvec,nyvec,nzvec]), max([kx,ky,kz]));
val = zeros(nxvec,nyvec,nzvec);

for i = 1:nx+kx-1
    if (xknot(i) > xknot(i+1))
        fprintf('subroutine dbs3gd:\n')
        fprintf('xknot(i) <= xknot(i+1) required.\n')
        fprintf(' i=%i, xknot(i)=%f, xknot(i+1)=%f\n',i, xknot(i), xknot(i+1))
        fprintf('\n')
        fprintf(' xknot=%f\n',xknot)
        error('here')
    end
end

for i = 1:nxvec
    if ((xvec(i)<xknot(1)) || (xvec(i)>xknot(nx+kx)))
        fprintf('subroutine dbs3gd:\n')
        fprintf('ix with xknot(ix) <= x < xknot(ix+1) required.\n')
        fprintf('x = %f\n', xvec(i))
        error('here')
    end
end

leftx(1) = 0;

leftx(1) = huntn(xknot,nx+kx,kx,xvec(1),leftx(1));

for ix = 2:nxvec
    leftx(ix) = leftx(ix-1);
    same = (xknot(leftx(ix)) <= xvec(ix)) &&  (xvec(ix) <= xknot(leftx(ix)+1));
    if(~same )
        leftx(ix) = leftx(ix) + 1;
        next      = (xknot(leftx(ix)) <= xvec(ix)) &&  (xvec(ix) <= xknot(leftx(ix)+1));
        if (~next)
            leftx(ix) = huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix));
        end
    end
end

for i = 1:ny+ky-1
    if (yknot(i) > yknot(i+1))
        fprintf('subroutine dbs3gd:\n')
        fprintf('yknot(i) <= yknot(i+1) required.\n')
        fprintf(' i=%i, yknot(i)=%f, yknot(i+1)=%f\n',i, yknot(i), yknot(i+1))
        fprintf('\n')
        fprintf(' yknot=%f\n',yknot)
        error('here')
    end
end

for i = 1:nyvec
    if ((yvec(i)<yknot(1)) || (yvec(i)>yknot(ny+ky)))
        fprintf('subroutine dbs3gd:\n')
        fprintf('iy with yknot(iy) <= y < yknot(iy+1) required.\n')
        fprintf('y = %f\n', yvec(i))
        error('here')
    end
end

lefty(1) = 0;

lefty(1) = huntn(yknot,ny+ky,ky,yvec(1),lefty(1));

for iy = 2:nyvec
    lefty(iy) = lefty(iy-1);
    same = (yknot(lefty(iy)) <= yvec(iy)) &&  (yvec(iy) <= yknot(lefty(iy)+1));
    if(~same )
        lefty(iy) = lefty(iy) + 1;
        next      = (yknot(lefty(iy)) <= yvec(iy)) &&  (yvec(iy) <= yknot(lefty(iy)+1));
        if (~next)
            lefty(iy) = huntn(yknot,ny+ky,ky,yvec(iy),lefty(iy));
        end
    end
end

for i = 1:nz+kz-1
    if (zknot(i) > zknot(i+1))
        fprintf('subroutine dbs3gd:\n')
        fprintf('zknot(i) <= zknot(i+1) required.\n')
        fprintf(' i=%i, zknot(i)=%f, zknot(i+1)=%f\n',i, zknot(i), zknot(i+1))
        fprintf('\n')
        fprintf(' zknot=%f\n',zknot)
        error('here')
    end
end

for i = 1:nzvec
    if ((zvec(i)<zknot(1)) || (zvec(i)>zknot(nz+kz)))
        fprintf('subroutine dbs3gd:\n')
        fprintf('iz with zknot(iz) <= z < zknot(iz+1) required.\n')
        fprintf('z = %f\n', zvec(i))
        error('here')
    end
end

leftz(1) = 0;

leftz(1) = huntn(zknot,nz+kz,kz,zvec(1),leftz(1));

for iz = 2:nzvec
    leftz(iz) = leftz(iz-1);
    same = (zknot(leftz(iz)) <= zvec(iz)) &&  (zvec(iz) <= zknot(leftz(iz)+1));
    if(~same )
        leftz(iz) = leftz(iz) + 1;
        next      = (zknot(leftz(iz)) <= zvec(iz)) &&  (zvec(iz) <= zknot(leftz(iz)+1));
        if (~next)
            leftz(iz) = huntn(zknot,nz+kz,kz,zvec(iz),leftz(iz));
        end
    end
end

if ((iderx == 0)  &&  (idery == 0)  &&  (iderz ==0))
    
    for ix = 1:nxvec
        biatx(ix,1) = 1.0;
    end
    
    for ik = 1:kx-1
        for ix = 1:nxvec
            dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix);
            dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik);
            save1(ix) = 0;
        end
        
        for il = 1:ik
            for ix = 1:nxvec
                term(ix)     = biatx(ix,il) / (dr(ix,il) + dl(ix,ik+1-il));
                biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix);
                save1(ix)    = dl(ix,ik+1-il) * term(ix);
            end
        end
        
        for ix = 1:nxvec
            biatx(ix,ik+1) = save1(ix);
        end
    end
    
    for iy = 1:nyvec
        biaty(iy,1) = 1.0;
    end
    
    for ik = 1:ky-1
        for iy = 1:nyvec
            dr(iy,ik) = yknot(lefty(iy)+ik) - yvec(iy);
            dl(iy,ik) = yvec(iy) - yknot(lefty(iy)+1-ik);
            save1(iy) = 0;
        end
        
        for il = 1:ik
            for iy = 1:nyvec
                term(iy)     = biaty(iy,il) / (dr(iy,il) + dl(iy,ik+1-il));
                biaty(iy,il) = save1(iy) + dr(iy,il) * term(iy);
                save1(iy)    = dl(iy,ik+1-il) * term(iy);
            end
        end
        
        for iy = 1:nyvec
            biaty(iy,ik+1) = save1(iy);
        end
    end
    
    for iz = 1:nzvec
        biatz(iz,1) = 1;
    end
    
    for ik = 1:kz-1
        for iz = 1:nzvec
            dr(iz,ik) = zknot(leftz(iz)+ik) - zvec(iz);
            dl(iz,ik) = zvec(iz) - zknot(leftz(iz)+1-ik);
            save1(iz) = 0;
        end
        
        for il = 1:ik
            for iz = 1:nzvec
                term(iz)     = biatz(iz,il) / (dr(iz,il) + dl(iz,ik+1-il));
                biatz(iz,il) = save1(iz) + dr(iz,il) * term(iz);
                save1(iz)    = dl(iz,ik+1-il) * term(iz);
            end
        end
        
        for iz = 1:nzvec
            biatz(iz,ik+1) = save1(iz);
        end
    end
    
    for iz = 1:nzvec
        for iy = 1:nyvec
            for ix = 1:nxvec
                val(ix,iy,iz) = 0;
            end
        end
    end
    
    for ikz = 1:kz
        for iky = 1:ky
            for ikx = 1:kx
                for iz = 1:nzvec
                    for iy = 1:nyvec
                        for ix = 1:nxvec
                            val(ix,iy,iz) = val(ix,iy,iz) + biatx(ix,ikx) * biaty(iy,iky) * biatz(iz,ikz) * bcoef(leftx(ix)-kx+ikx,lefty(iy)-ky+iky,leftz(iz)-kz+ikz);
                        end
                    end
                end
            end
        end
    end
    
else
    
    for iz = 1:nzvec
        for iy = 1:nyvec
            for ix = 1:nxvec
                val(ix,iy,iz) = dbs3dr(iderx,idery,iderz,xvec(ix),yvec(iy),zvec(iz),kx,ky,kz,xknot,yknot,zknot,nx,ny,nz,bcoef);
            end
        end
    end
    
end
