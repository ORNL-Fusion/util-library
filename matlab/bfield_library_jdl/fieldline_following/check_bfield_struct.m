function ierr = check_bfield_struct(bfield)
% checks the fields in bfield structure
ierr = 0;
errstr = 'Error from check_bfield_struct: ';
if ~isstruct(bfield)
    fprintf([errstr,'bfield must be a structure.\n'])
    ierr = 1;
    return;
end
if ~isfield(bfield,'type')
    fprintf([errstr,'bfield structure must contain field "type"\n'])
    ierr = 1;
    return;
end
if ~ischar(bfield.type)
    fprintf([errstr,bfield.type must be a character array\n'])
    ierr = 1;
    return;
end
if ~isfield(bfield,'nsym')
    fprintf([errstr,'nsym required for all bfield structures\n'])
    ierr = 1;
    return;
end

switch bfield.type
    case 'gfile'
        if ~isfield(bfield,'g')
            fprintf([errstr,'for gfile type structure must contain field "g"\n'])
            ierr = 1;
            return;
        end
    case 'gfile+coils'
        errstr2 = 'for gfile+coils ';
        if ~isfield(bfield,'g')
            fprintf([errstr,errstr2,'type structure must contain field "g"\n'])
            ierr = 1;
            return;
        end        
        if ~isfield(bfield,'coil')
            fprintf([errstr,errstr2,'type structure must contain field "coil"\n'])
            ierr = 1;
            return;
        end                
        if ~isfield(bfield,'current')
            fprintf([errstr,errstr2,'type structure must contain field "current"\n'])
            ierr = 1;
            return;
        end  
    case 'vmec'
        errstr2 = 'for vmec ';
        if ~isfield(bfield,'wout')
            fprintf([errstr,errstr2,'type structure must contain field "wout"\n'])
            ierr = 1;
            return;
        end        
    case 'just_coils'
        errstr2 = 'for just_coils ';    
        if ~isfield(bfield,'coil')
            fprintf([errstr,errstr2,'type structure must contain field "coil"\n'])
            ierr = 1;
            return;
        end                
        if ~isfield(bfield,'current')
            fprintf([errstr,errstr2,'type structure must contain field "current"\n'])
            ierr = 1;
            return;
        end      
    case 'MPEX'
        errstr2 = 'for MPEX ';    
        if ~isfield(bfield,'coil')
            fprintf([errstr,errstr2,'type structure must contain field "coil"\n'])
            ierr = 1;
            return;
        end                
        if ~isfield(bfield,'current')
            fprintf([errstr,errstr2,'type structure must contain field "current"\n'])
            ierr = 1;
            return;
        end            
    case 'ipec_eq'
        errstr2 = 'for ipec_eq ';
        if ~isfield(bfield,'ipec')
            fprintf([errstr,errstr2,'type structure must contain field "ipec"\n'])
            ierr = 1;
            return;
        end         
        if ~isfield(bfield.ipec,'eq')
            fprintf([errstr,errstr2,'type structure must contain field "ipec.eq"\n'])
            ierr = 1;
            return;
        end        
    case {'ipec_vac','ipec_vac_only'}
        errstr2 = 'for ipec_vac ';
        if ~isfield(bfield,'ipec')
            fprintf([errstr,errstr2,'type structure must contain field "ipec"\n'])
            ierr = 1;
            return;
        end         
        if ~isfield(bfield.ipec,'vac')
            fprintf([errstr,errstr2,'type structure must contain field "ipec.vac"\n'])
            ierr = 1;
            return;
        end          
    case {'ipec_pert','ipec_pert_only'}
        errstr2 = 'for ipec_pert ';
        if ~isfield(bfield,'ipec')
            fprintf([errstr,errstr2,'type structure must contain field "ipec"\n'])
            ierr = 1;
            return;
        end         
        if ~isfield(bfield.ipec,'pert')
            fprintf([errstr,errstr2,'type structure must contain field "ipec.pert"\n'])
            ierr = 1;
            return;
        end        
    case {'xpand_pert'}
        errstr2 = 'for xpand_pert ';
        if ~isfield(bfield,'xpand')
            fprintf([errstr,errstr2,'type structure must contain field "xpand"\n'])
            ierr = 1;
            return;
        end         
        if ~isfield(bfield.xpand,'Br')
            fprintf([errstr,errstr2,'type structure must contain field "xpand.Br"\n'])
            ierr = 1;
            return;
        end            
    case {'xpand_vac'}
        errstr2 = 'for xpand_vac ';
        if ~isfield(bfield,'xpand')
            fprintf([errstr,errstr2,'type structure must contain field "g"\n'])
            ierr = 1;
            return;
        end                
        if ~isfield(bfield,'xpand')
            fprintf([errstr,errstr2,'type structure must contain field "xpand"\n'])
            ierr = 1;
            return;
        end         
        if ~isfield(bfield.xpand,'Br')
            fprintf([errstr,errstr2,'type structure must contain field "xpand.Br"\n'])
            ierr = 1;
            return;
        end            
    case {'Bgrid'}
        errstr2 = 'for Bgrid';
        if ~isfield(bfield,'Bgrid')
            fprintf([errstr,errstr2,'type structure must contain field "Bgrid"\n'])
            ierr = 1;
            return;
        end  
    case {'Aspline'}
        errstr2 = 'for Aspline';
        if ~isfield(bfield,'Arcoeff')
            fprintf([errstr,errstr2,'type structure must contain field "Arcoeff"\n'])
            ierr = 1;
            return;
        end            
    otherwise
        fprintf([errstr,'Did not recognize bfield type\n'])
        fprintf('Supported types are:\n')
        fprintf('     gfile\n')
        fprintf('     gfile+coils\n')
        fprintf('     vmec\n')        
        fprintf('     just_coils\n')        
        fprintf('     MPEX\n')      
        fprintf('     ipec_eq\n')     
        fprintf('     ipec_vac\n')     
        fprintf('     ipec_pert\n')    
        fprintf('     ipec_vac_only\n')     
        fprintf('     ipec_pert_only\n')        
        fprintf('     xpand_pert\n') 
        fprintf('     xpand_vac\n') 
        fprintf('     Bgrid\n')
        ierr = 1;
        error('Throwing error')        
        return;
end
