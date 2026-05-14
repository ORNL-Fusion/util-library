function [CoilGeometry,currentPerWinding] = setup_MPEX_coils(current_in,verbose)
% [CoilGeometry,currentPerWinding] = setup_MPEX_coils(current_in,verbose)
% Sets the arrays cur for each winding of the coil set.
%
% current_in can be specified in two ways
%
% If current_in is a string or character vector matching one of the
% configuration names in get_MPEX_currents_by_config, the corresponding
% coil-current column will be loaded.
%
% If current_in is a numeric array with length equal to the number of coils
% in define_MPEX_coil_filaments (23), then the values will be used as the
% winding current for each coil [Amps].
%
% verbose ~= 0 will print the coil currents to screen.

% current_in = the current_in in each WINDING (amps)
% J.D. Lore
CoilGeometry = define_MPEX_coils;

% Current configuration data is supplied by the private MPEX-modeling-data repo.
[config_name,data] = get_MPEX_currents_by_config;


if nargin == 0
    error('Must specify inputs')
end
if nargin < 2
    verbose = 0;
end

currentPerWinding = zeros(1,CoilGeometry.ncoils); %#ok<PREALL>
if ischar(current_in) || (isstring(current_in) && isscalar(current_in))
    current_name = char(string(current_in));
    ind = find(strcmp(current_name,config_name),1,'first');
    if isempty(ind)
        error('Unknown MPEX current configuration "%s"',current_name)
    end
    currentPerWinding = data(:,ind).';
elseif isnumeric(current_in) && isvector(current_in) && length(current_in) == CoilGeometry.ncoils
    currentPerWinding = current_in;
else
    error('current_in must be a config name or a numeric array of length fil.ncoils = %d',CoilGeometry.ncoils)
end

if verbose
    fprintf('------------------------------------------------------------------------------------------------\n')
    fprintf('Coil currents set to: \n')
    fprintf('%8d',1:CoilGeometry.ncoils)
    fprintf('\n')
    fprintf('%8.1f',currentPerWinding)
    fprintf('\n')
    if ischar(current_in) || (isstring(current_in) && isscalar(current_in))
        fprintf('Using configuration %s\n',current_name)
    end
    fprintf('------------------------------------------------------------------------------------------------\n')
end
