function [dydx,ierr] = choose_fl_derivs_dphi(x,y,bfield,nowarn)
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
    case 'ipec_eq'        
        [dydx,ierr] = fl_derivs_dphi_ipec(x,y,bfield,nowarn,1);        
    case 'ipec_vac'        
        [dydx,ierr] = fl_derivs_dphi_ipec(x,y,bfield,nowarn,2);                
    otherwise
        fprintf('Did not recognize bfield type\n')
        fprintf('Supported types are:\n')
        fprintf('     gfile\n')
        fprintf('     gfile+coils\n')
        fprintf('     vmec\n')        
        fprintf('     just_coils\n')
        fprintf('     ipec_eq\n')
        fprintf('     ipec_vac\n')
        ierr = 1;
        return;
end