function rmp = build_d3d_ccoils_jl(taper,ntorpts)
%
%  SET UP SO TAPER SHOULD BE TAKEN DIRECTLY FROM curntc in diiidsup.in
%
% Builds the coil geometry and sets the current (in Amp-turns) for 
% the 6 D3D C coils. A single loop is returned for each entire 
% coil set, with one turn per coil.
%
% Coils 1:6 --> C79, C139, C199, C259, C319, C19
%   Coil toroidal angles are in RH cylindrical coordinates, so negative
%   of machine geographical toroidal angle.  <-- Following M. Shaffer routines
%
%
% Inputs:
%	taper(1:6) -- Contains the coil current in Amps for the coils. This
%                 value is multiplied by number of turns/coil.
%	ntorpts    -- Number of toroidal cuts used to define each coil.  (default 4)
%
% Output:
%	rmp.coil(npts,1:3) -- Coil set coordinates in meters. The second index corresponds
%	                      to the x,y, or z coordinate in meters. 
%   rmp.current(npts)  -- Current in each filament. The filaments connecting the RWM 
%					      coils have 0 current, and the last value is irrelevant 
%					     (number of filaments = number of points - 1).
%
% JDL

if nargin < 2 
    ntorpts = 4;
end
npts = 2*ntorpts + 1;

% Coils are approximately 'picture frame' centered at midplane
num_coils = 6;
phicens = [-79,-139,-199,-259,-319,-19]*pi/180;  % Center of each coil in degrees
phiext = 58*pi/180;						  % Toroidal extent of each coil in degrees
R = [3.23,3.23];                          % Major radius of upper and lower toroidal coil arc (meters)
Z = [0.8,-0.8];                     % Vertical position of coil (meters)
nturn = 4;  % number of electrical turns

coil = zeros(num_coils*npts,3);
current = zeros(num_coils*npts,1);
taper = nturn*taper;

for i = 0:(num_coils-1)
    phi = phiext*[0:ntorpts-1]/(ntorpts-1) + phicens(i+1) - phiext/2;

    coil(i*npts+1:i*npts+ntorpts,1) = R(1)*cos(phi);
    coil(i*npts+1:i*npts+ntorpts,2) = R(1)*sin(phi);
    coil(i*npts+1:i*npts+ntorpts,3) = Z(1);
    coil(i*npts+ntorpts+1:i*npts+2*ntorpts,1) = R(2)*cos(phi(end:-1:1));
    coil(i*npts+ntorpts+1:i*npts+2*ntorpts,2) = R(2)*sin(phi(end:-1:1));
    coil(i*npts+ntorpts+1:i*npts+2*ntorpts,3) = Z(2);
    coil((i+1)*npts,:) = coil(i*npts+1,:);
    current(i*npts+1:(i+1)*npts) = taper(i+1);
    current((i+1)*npts) = 0.0;  % Sticks connecting the coils have no current.
end
rmp.coil = coil;
rmp.current = current;

