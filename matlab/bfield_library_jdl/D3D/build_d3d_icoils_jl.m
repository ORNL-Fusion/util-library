function rmp = build_d3d_icoils_jl(taper,ntorpts)
%
%  SET UP SO TAPER SHOULD BE TAKEN DIRECTLY FROM curntic in diiidsup.in
%
% Builds the coil geometry and sets the current (in Amp-turns) for 
% the 12 D3D I coils. A single loop is returned for each entire 
% coil set, with one turn per coil.
%
% Coils 1:12 -->  IU30, IU90, IU150, IU210, IU270, IU330, IL30, IL90, IL150, IL210, IL270, IL330.
%   Coil toroidal angles are in RH cylindrical coordinates, so negative
%   of machine geographical toroidal angle.  <-- Following M. Shaffer routines
%
%
% Inputs:
%	taper(1:12) -- Contains the coil current in Amps for the coils. This
%                 value is multiplied by number of turns/coil.
%	ntorpts    -- Number of toroidal cuts used to define each coil.  (default 6)
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
    ntorpts = 6;
end
npts = 2*ntorpts + 1;

% Coils are approximately 'picture frame' 
num_coils = 12;
phicens = [-32.7,-87.3,-152.7,-207.3,-272.7,-327.3,...
           -32.7,-87.3,-152.7,-207.3,-272.7,-327.3]*pi/180;  % Center of each coil in degrees
phiext = 51.72*pi/180;      % Toroidal extent of each coil
% 2003 
R = [2.184,2.394];          % Major radius of upper and lower toroidal coil arc (meters)
Z = [1.012,0.504];          % Vertical position of coil (meters)  --- Upper coil, lower just invert and flip Z
% % 2006 'revised'
% R = [2.164,2.373];          % Major radius of upper and lower toroidal coil arc (meters)
% Z = [1.016,0.504];          % Vertical position of coil (meters)  --- Upper coil, lower just invert and flip Z
nturn = 1;                  % number of electrical turns


coil = zeros(num_coils*npts,3);
current = zeros(num_coils*npts,1);
taper = -nturn*taper;    % The - sign here accounts for DIII-D convention!

% UPPER COILS
for i = 0:num_coils/2 -1
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
% LOWER COILS
for i = num_coils/2:(num_coils - 1)
    phi = phiext*[0:ntorpts-1]/(ntorpts-1) + phicens(i+1) - phiext/2;
    coil(i*npts+1:i*npts+ntorpts,1) = R(2)*cos(phi);
    coil(i*npts+1:i*npts+ntorpts,2) = R(2)*sin(phi);
    coil(i*npts+1:i*npts+ntorpts,3) = -Z(2);
    coil(i*npts+ntorpts+1:i*npts+2*ntorpts,1) = R(1)*cos(phi(end:-1:1));
    coil(i*npts+ntorpts+1:i*npts+2*ntorpts,2) = R(1)*sin(phi(end:-1:1));
    coil(i*npts+ntorpts+1:i*npts+2*ntorpts,3) = -Z(1);
    coil((i+1)*npts,:) = coil(i*npts+1,:);
    current(i*npts+1:(i+1)*npts) = taper(i+1);
    current((i+1)*npts) = 0.0;  % Sticks connecting the coils have no current.
end
rmp.coil = coil;
rmp.current = current;

