function [dydx,ierr] = choose_fl_derivs_dl(y,bfield,nowarn)
if nargin < 3
    nowarn = 0;
end

switch bfield.type
    case 'gfile'
        [dydx,ierr] = fl_derivs_dl_gfile(y,bfield,nowarn);
    case 'gfile+coils'
        [dydx,ierr] = fl_derivs_dl_gfile_coils(y,bfield,nowarn);
    case 'just_coils'        
        [dydx,ierr] = fl_derivs_dl_just_coils(y,bfield,nowarn);     
    case 'MPEX'
        [dydx,ierr] = fl_derivs_dl_MPEX(y,bfield,nowarn);     
    otherwise
        fprintf([errstr,'Did not recognize bfield type\n'])
        fprintf('Supported types are:\n')
        fprintf('     gfile\n')
        fprintf('     gfile+coils\n')
        fprintf('     just_coils\n')        
        fprintf('     MPEX\n')        
        ierr = 1;
        return;
end