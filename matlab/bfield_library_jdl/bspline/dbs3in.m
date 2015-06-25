function bcoef=dbs3in(nx,xvec,ny,yvec,nz,zvec,xyzdata,ldf,mdf,kx,ky,kz,xknot,yknot,zknot)

% !
% !  Computes a three-dimensional tensor-product spline interpolant,
% !  returning the tensor-product B-spline coefficients.
% !
% !   nx      - number of data points in the x-direction.  (input)
% !   xvec    - array of length nxdata containing the data points in
% !             the x-direction.  (input)
% !             xdata must be increasing.
% !   ny      - number of data points in the y-direction.  (input)
% !   yvec    - array of length nydata containing the data points in
% !             the y-direction.  (input)
% !             ydata must be increasing.
% !   nz      - number of data points in the z-direction.  (input)
% !   zvec    - array of length nzdata containing the data points in
% !             the z-direction.  (input)
% !             zdata must be increasing.
% !   xyzdata - array of size nx by ny by nz containing the
% !             values to be interpolated.  (input)
% !             xyzdata(i,j,k) contains the value at
% !             (xvec(i),yvec(j),zvec(k)).
% !   ldf     - leading dimension of fdata exactly as specified in the
% !             dimension statement of the calling program.  (input)
% !   mdf     - middle dimension of fdata exactly as specified in the
% !             dimension statement of the calling program.  (input)
% !   kx      - order of the spline in the x-direction.  (input)
% !             kxord must be less than or equal to nxdata.
% !   ky      - order of the spline in the y-direction.  (input)
% !             kyord must be less than or equal to nydata.
% !   kz      - order of the spline in the z-direction.  (input)
% !             kzord must be less than or equal to nzdata.
% !   xknot   - array of length nx+kx containing the knot
% !             sequence in the x-direction.  (input)
% !             xknot must be nondecreasing.
% !   yknot   - array of length ny+ky containing the knot
% !             sequence in the y-direction.  (input)
% !             yknot must be nondecreasing.
% !   zknot   - array of length nz+kz containing the knot
% !             sequence in the z-direction.  (input)
% !             zknot must be nondecreasing.
% !   bcoef   - array of length nx*ny*nz containing the
% !             tensor-product B-spline coefficients.  (output)
% !             bscoef is treated internally as a matrix of size nx
% !             by ny by nz.
% !

%     integer, intent(in) :: nx, ny, nz, kx, ky, kz
%     integer, intent(in) :: ldf, mdf
% 
%     real(kind=dbl), dimension(nx), intent(in)         :: xvec
%     real(kind=dbl), dimension(ny), intent(in)         :: yvec
%     real(kind=dbl), dimension(nz), intent(in)         :: zvec
%     real(kind=dbl), dimension(nx+kx), intent(in)      :: xknot
%     real(kind=dbl), dimension(ny+ky), intent(in)      :: yknot
%     real(kind=dbl), dimension(nz+kz), intent(in)      :: zknot
%     real(kind=dbl), dimension(ldf,mdf,nz), intent(in) :: xyzdata
%     real(kind=dbl), dimension(nx,ny,nz), intent(out)  :: bcoef
% 
%     integer                                :: iz
%     !real(kind=dbl), dimension(nx,ny,nz)    :: work1
%     real(kind=dbl), dimension(nz)          :: work2
%     real(kind=dbl), dimension((2*kz-1)*nz) :: work3
%     real(kind=dbl), dimension(:,:,:), allocatable :: work1


%     allocate (work1(nx,ny,nz))
    [~,~,work1] = spli3d(zvec,ldf,mdf,xyzdata,zknot,nz,kz,nx,ny,nx,ny,nz);
    bcoef = zeros(nx,ny,nz);
    for iz = 1:nz
       bcoef(:,:,iz) = dbs2in(nx,xvec,ny,yvec,work1(:,:,iz),nx,kx,ky,xknot,yknot);
    end
