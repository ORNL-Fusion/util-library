function [Bout,ierr] = bfield_general_rzphi(R,Z,phi_radian,bfield,nowarn)
% [Bout,ierr] = bfield_general_rzphi(R,Z,phi_radian,bfield,nowarn)
if nargin < 5
    nowarn = 0;
end

switch bfield.type
    case 'gfile'        
        [Bout,ierr] = bfield_geq_bicub(bfield.g,R,Z,nowarn);
    case 'gfile+coils'
        [Bout,ierr] = bfield_geq_bicub(bfield.g,R,Z,nowarn);
        [Br,Bphi,Bz]=bfield_bs_cyl(R,phi_radian,Z,bfield.coil,bfield.current,nowarn);
        Bout.br = Bout.br + Br;
        Bout.bphi = Bout.bphi + Bphi;
        Bout.bz = Bout.bz + Bz;
%     case 'vmec'        
%         [dydx,ierr] = fl_derivs_dphi_vmec(x,y,bfield,nowarn);
    case 'just_coils'        
        [Br,Bphi,Bz]=bfield_bs_cyl(R,phi_radian,Z,bfield.coil,bfield.current,nowarn);
        Bout.br = Br;
        Bout.bphi = Bphi;
        Bout.bz = Bz;
        ierr = 0;
%         [dydx,ierr] = fl_derivs_dphi_just_coils(x,y,bfield,nowarn);
%     case 'ipec_eq'        
%         [dydx,ierr] = fl_derivs_dphi_ipec(x,y,bfield,nowarn,0);        
%     case 'ipec_vac'        
%         [dydx,ierr] = fl_derivs_dphi_ipec(x,y,bfield,nowarn,1);                
%     case 'ipec_pert'        
%         [dydx,ierr] = fl_derivs_dphi_ipec(x,y,bfield,nowarn,2);        
%     case 'ipec_vac_only'        
%         [dydx,ierr] = fl_derivs_dphi_ipec(x,y,bfield,nowarn,3);                
%     case 'ipec_pert_only'        
%         [dydx,ierr] = fl_derivs_dphi_ipec(x,y,bfield,nowarn,4);               
%     case 'xpand_pert'        
%         [dydx,ierr] = fl_derivs_dphi_xpand(x,y,bfield,nowarn,0);                
%     case 'xpand_vac'        
%         [dydx,ierr] = fl_derivs_dphi_xpand(x,y,bfield,nowarn,1);      
    case 'Bgrid'        
        [Br,Bz,Bphi]=bfield_grid(R,Z,phi_radian,bfield.Bgrid,nowarn);        
        Bout.br = Br;
        Bout.bphi = Bphi;
        Bout.bz = Bz;
        ierr = 0;
    case 'Aspline'
        [Br,Bz,Bphi,ierr]=bfield_Aspline(R,Z,phi_radian,bfield.Arcoeff,bfield.Azcoeff,bfield.Aphicoeff,bfield.spline_info);              
        Bout.br = Br;
        Bout.bphi = Bphi;
        Bout.bz = Bz;
    otherwise
        fprintf('Did not recognize bfield type\n')
        fprintf('Supported types are:\n')
        fprintf('     gfile\n')
        fprintf('     gfile+coils\n')
%         fprintf('     vmec\n')        
        fprintf('     just_coils\n')
%         fprintf('     ipec_eq\n')
%         fprintf('     ipec_vac\n')
%         fprintf('     ipec_pert\n')
%         fprintf('     ipec_vac_only\n')
%         fprintf('     ipec_pert_only\n')        
%         fprintf('     xpand_vac\n')
%         fprintf('     xpand_pert\n')
        fprintf('     Bgrid\n')
        fprintf('     Agrid\n')
        error('quitting')
        ierr = 1;
        return;
end