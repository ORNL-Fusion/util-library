function dbs3ine
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % IMSL name:  dbs3in (double precision version)
% %
% % purpose:    compute a three-dimensional tensor-product spline
% %             interpolant, returning the tensor-product B-spline
% %             coefficients.
% %
% % usage:      call dbs3in(nxdata, xdata, nydata, ydata, nzdata,
% %                         zdata, fdata, ldf, mdf, kxord, kyord,
% %                         kzord, xknot, yknot, zknot, bscoef)
% %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %
% %        specifications for parameters
% %
% 
%       integer    kxord, kyord, kzord, ldf, mdf, nxdata, nxknot, nxvec,  &
%      &           nydata, nyknot, nyvec, nzdata, nzknot, nzvec
%       parameter  (kxord=5, kyord=2, kzord=3, nxdata=21, nxvec=4,        &
%      &           nydata=6, nyvec=4, nzdata=8, nzvec=2, ldf=nxdata,      &
%      &           mdf=nydata, nxknot=nxdata+kxord, nyknot=nydata+kyord,  &
%      &           nzknot=nzdata+kzord)
% 
%       integer    i, j, k, nxcoef, nycoef, nzcoef
%       double precision  bscoef(nxdata,nydata,nzdata), f,                &
%      &           fdata(ldf,mdf,nzdata), value(nxvec,nyvec,nzvec)        &
%      &           , x, xdata(nxdata), xknot(nxknot), xvec(nxvec), y,     &
%      &           ydata(nydata), yknot(nyknot), yvec(nyvec), z,          &
%      &           zdata(nzdata), zknot(nzknot), zvec(nzvec)
% 
% %
% %        define function.
% %

kxord = 5;
kyord=2;
kzord=3;
nxdata=21;
nxvec=4;
nydata=6; 
nyvec=4; 
nzdata=8; 
nzvec=2; 
ldf=nxdata;
mdf=nydata; 
nxknot=nxdata+kxord; 
nyknot=nydata+kyord;
nzknot=nzdata+kzord;


% %
% %        set up x-interpolation points
% %

      for i = 1:nxdata
         xdata(i) = (i-11)/10.0;
      end

%
%        set up y-interpolation points
%

      for i = 1:nydata
         ydata(i) = (i-1)/(nydata-1);
      end

%
%        set up z-interpolation points
%

      for i = 1:nzdata
         zdata(i) = (i-1)/(nzdata-1);
      end

%
%        generate knots
%

      xknot= dbsnak (nxdata, xdata, kxord);
      yknot = dbsnak (nydata, ydata, kyord);
      zknot = dbsnak (nzdata, zdata, kzord);

%
%        generate fdata
%

      for k = 1:nzdata
         for i = 1:nydata
            for j = 1:nxdata
               fdata(j,i,k) = f(xdata(j),ydata(i),zdata(k));
            end
         end
      end

%
%        interpolate
%

      bscoef = dbs3in (nxdata, xdata, nydata, ydata, nzdata, zdata, fdata, ldf, mdf, kxord, kyord, kzord, xknot, yknot, zknot);

      nxcoef = nxdata;
      nycoef = nydata;
      nzcoef = nzdata;

%
%        write heading
%

fprintf('   x        y     z    s(x,y,z)    error\n')

%        print over a grid of [-1.0,1.0] x [0.0,1.0] x [0.0,1.0]
%        at 32 points.

      for i = 1:nxvec
         xvec(i) = 2.0*((i-1)/3.0) - 1.0;
      end

      for i = 1:nyvec
         yvec(i) = (i-1)/3.0;
      end

      for i = 1:nzvec
         zvec(i) = (i-1);
      end

%
%        call the evaluation routine.
%

      value = dbs3gd (0, 0, 0, nxvec, xvec, nyvec, yvec, nzvec, zvec, kxord, kyord, kzord, xknot, yknot, zknot, nxcoef, nycoef, nzcoef, bscoef, nxvec, nyvec);

      for i = 1:nxvec
         for j = 1:nyvec
            for k = 1:nzvec
               fprintf('%13.4f %13.4f %13.4f %13.4f %13.6f\n', xvec(i), yvec(k), zvec(k), value(i,j,k), f(xvec(i),yvec(j),zvec(k)) - value(i,j,k));
            end
         end
      end

     function f = f(x,y,z)
           f = x*x*x + x*y*z;
     end

      end
% 
% %           x           y          z          s(x,y,z)       error
% %       -1.0000       0.0000       0.0000      -1.0000     0.000000
% %       -1.0000       0.3333       1.0000      -1.0000     0.000000
% %       -1.0000       0.0000       0.0000      -1.0000     0.000000
% %       -1.0000       0.3333       1.0000      -1.3333     0.000000
% %       -1.0000       0.0000       0.0000      -1.0000     0.000000
% %       -1.0000       0.3333       1.0000      -1.6667     0.000000
% %       -1.0000       0.0000       0.0000      -1.0000     0.000000
% %       -1.0000       0.3333       1.0000      -2.0000     0.000000
% %       -0.3333       0.0000       0.0000      -0.0370     0.000000
% %       -0.3333       0.3333       1.0000      -0.0370     0.000000
% %       -0.3333       0.0000       0.0000      -0.0370     0.000000
% %       -0.3333       0.3333       1.0000      -0.1481     0.000000
% %       -0.3333       0.0000       0.0000      -0.0370     0.000000
% %       -0.3333       0.3333       1.0000      -0.2593     0.000000
% %       -0.3333       0.0000       0.0000      -0.0370     0.000000
% %       -0.3333       0.3333       1.0000      -0.3704     0.000000
% %        0.3333       0.0000       0.0000       0.0370     0.000000
% %        0.3333       0.3333       1.0000       0.0370     0.000000
% %        0.3333       0.0000       0.0000       0.0370     0.000000
% %        0.3333       0.3333       1.0000       0.1481     0.000000
% %        0.3333       0.0000       0.0000       0.0370     0.000000
% %        0.3333       0.3333       1.0000       0.2593     0.000000
% %        0.3333       0.0000       0.0000       0.0370     0.000000
% %        0.3333       0.3333       1.0000       0.3704     0.000000
% %        1.0000       0.0000       0.0000       1.0000     0.000000
% %        1.0000       0.3333       1.0000       1.0000     0.000000
% %        1.0000       0.0000       0.0000       1.0000     0.000000
% %        1.0000       0.3333       1.0000       1.3333     0.000000
% %        1.0000       0.0000       0.0000       1.0000     0.000000
% %        1.0000       0.3333       1.0000       1.6667     0.000000
% %        1.0000       0.0000       0.0000       1.0000     0.000000
% %        1.0000       0.3333       1.0000       2.0000     0.000000
