function [dydx,ierr] = choose_fl_derivs(x,y,bfield,nowarn)
if nargin < 4
    nowarn = 0;
end

switch bfield.type
    case 'gfile'
        [dydx,ierr] = fl_derivs_dphi_gfile(y,bfield,nowarn);
    case 'gfile+coils'
        [dydx,ierr] = fl_derivs_dphi_gfile_coils(x,y,bfield,nowarn);
    case 'vmec'        
        [dydx,ierr] = fl_derivs_dphi_vmec(x,y,bfield,nowarn);
    case 'just_coils'        
        [dydx,ierr] = fl_derivs_dphi_just_coils(x,y,bfield,nowarn);
    otherwise
        fprintf([errstr,'Did not recognize bfield type\n'])
        fprintf('Supported types are:\n')
        fprintf('     gfile\n')
        fprintf('     gfile+coils\n')
        fprintf('     vmec\n')        
        fprintf('     just_coils\n')
        ierr = 1;
        return;
end