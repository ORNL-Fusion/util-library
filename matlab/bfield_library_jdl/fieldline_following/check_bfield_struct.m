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
            fprintf([errstr,errstr2,'type structure must contain field "g"\n'])
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
    otherwise
        fprintf([errstr,'Did not recognize bfield type\n'])
        fprintf('Supported types are:\n')
        fprintf('     gfile\n')
        fprintf('     gfile+coils\n')
        fprintf('     vmec\n')        
        fprintf('     just_coils\n')        
        fprintf('     MPEX\n')      
        error('Throwing error')
        ierr = 1;
        return;
end
