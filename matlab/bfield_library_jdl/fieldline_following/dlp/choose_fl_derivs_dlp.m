function [dydx,ierr] = choose_fl_derivs_dlp(y,bfield,nowarn)
if nargin < 3
    nowarn = 0;
end

switch bfield.type
    case 'gfile'
        [dydx,ierr] = fl_derivs_dlp_gfile(y,bfield,nowarn);
    otherwise
        fprintf('Did not recognize bfield type\n')
        fprintf('Supported types are:\n')
        fprintf('     gfile\n')
        error('Add more bfield types!')
end