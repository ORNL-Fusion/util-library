function [dydx,ierr] = choose_fl_derivs_dz(x,y,bfield,nowarn)
if nargin < 4
    nowarn = 0;
end

switch bfield.type
    case 'just_coils'        
        [dydx,ierr] = fl_derivs_dz_just_coils(x,y,bfield,nowarn);     
    case 'MPEX'
        [dydx,ierr] = fl_derivs_dz_MPEX(x,y,bfield,nowarn);     
    otherwise
        fprintf([errstr,'Did not recognize bfield type\n'])
        fprintf('Supported types are:\n')
        fprintf('     just_coils\n')        
        fprintf('     MPEX\n')        
        ierr = 1;
        return;
end