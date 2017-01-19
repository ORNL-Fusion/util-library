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
        [dydx,ierr] = fl_derivs_dphi_ipec(x,y,bfield,nowarn,0);        
    case 'ipec_vac'        
        [dydx,ierr] = fl_derivs_dphi_ipec(x,y,bfield,nowarn,1);                
    case 'ipec_pert'        
        [dydx,ierr] = fl_derivs_dphi_ipec(x,y,bfield,nowarn,2);        
    case 'ipec_vac_only'        
        [dydx,ierr] = fl_derivs_dphi_ipec(x,y,bfield,nowarn,3);                
    case 'ipec_pert_only'        
        [dydx,ierr] = fl_derivs_dphi_ipec(x,y,bfield,nowarn,4);               
    case 'xpand_pert'        
        [dydx,ierr] = fl_derivs_dphi_xpand(x,y,bfield,nowarn,0);                
    case 'xpand_vac'        
        [dydx,ierr] = fl_derivs_dphi_xpand(x,y,bfield,nowarn,1);                        
    case 'Bgrid'        
        [dydx,ierr] = fl_derivs_dphi_Bgrid(x,y,bfield,nowarn);               
    case 'Aspline'        
        [dydx,ierr] = fl_derivs_dphi_Aspline(x,y,bfield,nowarn);         
    otherwise
        fprintf('Did not recognize bfield type\n')
        fprintf('Supported types are:\n')
        fprintf('     gfile\n')
        fprintf('     gfile+coils\n')
        fprintf('     vmec\n')        
        fprintf('     just_coils\n')
        fprintf('     ipec_eq\n')
        fprintf('     ipec_vac\n')
        fprintf('     ipec_pert\n')
        fprintf('     ipec_vac_only\n')
        fprintf('     ipec_pert_only\n')        
        fprintf('     xpand_vac\n')
        fprintf('     xpand_pert\n')
        fprintf('     Bgrid\n')
        fprintf('     Aspline\n')
        ierr = 1;
        error('barfing')
        return;
end