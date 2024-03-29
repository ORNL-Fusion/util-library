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
    case 'just_coils'        
        [Br,Bphi,Bz]=bfield_bs_cyl(R,phi_radian,Z,bfield.coil,bfield.current,nowarn);
        Bout.br = Br;
        Bout.bphi = Bphi;
        Bout.bz = Bz;
        ierr = 0; 
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
    case 'Bspline'
        [Br,Bz,Bphi,ierr]=bfield_Bspline(R,Z,phi_radian,bfield.Brcoeff,bfield.Bzcoeff,bfield.Bphicoeff,bfield.spline_info);              
        Bout.br = Br;
        Bout.bphi = Bphi;
        Bout.bz = Bz;        
    case 'MPEX'
        [Br,Bz] = bfield_circular_coils(bfield.coil,bfield.current,R,Z);
        Bout.br = Br;
        Bout.bz = Bz;
        Bout.bphi = zeros(size(Br));
        ierr = 0;
    case 'xdr'
        [Br,Bz,Bphi,ierr] = bint_xdr(R,Z,phi_radian,bfield.xdr,nowarn);
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
        fprintf('     Aspline\n')
        fprintf('     Bspline\n')        
        fprintf('      MPEX\n')
        fprintf('      xdr\n')
        error('quitting')
        ierr = 1;
        return;
end


% Bx = Br.*cp - Bphi.*sp;
% By = Br.*sp + Bphi.*cp;
% Bpol = = sqrt(B.bz.^2 + B.br.^2)