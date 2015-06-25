function bcoef=dbs2in(nx,xvec,ny,yvec,xydata,ldf,kx,ky,xknot,yknot)

% !
% !  Computes a two-dimensional tensor-product spline interpolant,
% !  returning the tensor-product B-spline coefficients.
% !
% !    nx     - number of data points in the x-direction.  (input)
% !    xvec   - array of length nx containing the data points in
% !             the x-direction.  (input)
% !             xdata must be strictly increasing.
% !    ny     - number of data points in the y-direction.  (input)
% !    yvec   - array of length ny containing the data points in
% !             the y-direction.  (input)
% !             ydata must be strictly increasing.
% !    xydata - array of size nx by nydata containing the values to
% !             be interpolated.  (input)
% !             fdata(i,j) is the value at (xdata(i),ydata(j)).
% !    ldf    - the leading dimension of fdata exactly as specified in
% !             the dimension statement of the calling program.
% !             (input)
% !    kx     - order of the spline in the x-direction.  (input)
% !             kxord must be less than or equal to nxdata.
% !    ky     - order of the spline in the y-direction.  (input)
% !             kyord must be less than or equal to nydata.
% !    xknot  - array of length nx+kx containing the knot
% !             sequence in the x-direction.  (input)
% !             xknot must be nondecreasing.
% !    yknot  - array of length ny+ky containing the knot
% !             sequence in the y-direction.  (input)
% !             yknot must be nondecreasing.
% !    bcoef  - array of length nx*ny containing the
% !             tensor-product B-spline coefficients.  (output)
% !             bscoef is treated internally as a matrix of size nxdata
% !             by nydata.
% !

%     integer, intent(in)                           :: nx, ny, kx, ky, ldf
% 
%     real(kind=dbl), dimension(nx), intent(in)     :: xvec
%     real(kind=dbl), dimension(ny), intent(in)     :: yvec
%     real(kind=dbl), dimension(nx+kx), intent(in)  :: xknot
%     real(kind=dbl), dimension(ny+ky), intent(in)  :: yknot
%     real(kind=dbl), dimension(ldf,*), intent(in)  :: xydata
%     real(kind=dbl), dimension(nx,ny), intent(out) :: bcoef
% 
%     real(kind=dbl), dimension(max(nx,ny),max(nx,ny))        :: work1
%     real(kind=dbl), dimension(max(nx,ny))                   :: work2
%     real(kind=dbl), dimension(max((2*kx-1)*nx,(2*ky-1)*ny)) :: work3


    [~,~,work1] = spli2d(xvec,ldf,xydata,xknot,nx,kx,ny);
    [~,~,bcoef] = spli2d(yvec,ny, work1, yknot,ny,ky,nx);
