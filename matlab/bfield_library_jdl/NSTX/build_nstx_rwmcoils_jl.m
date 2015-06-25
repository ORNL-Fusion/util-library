function rmp = build_nstx_rwmcoils_jl(taper,ntorpts)
% Builds the coil geometry and sets the current (in Amp-turns) for 
% the 6 NSTX RWM coils. A single loop is returned for each entire 
% coil set, with one turn per RWM coil.
%
% Inputs:
%	taper(1:6) -- Contains the coil current in Amps for RWM#1-6. This
%                 value is multiplied by 2 turns/coil. If length(taper) == 1
%                 then this is the current for RWM1 and an n=3 phase is used.
%                 This should have the sign of the irwm# MDSplus signals.
%	ntorpts -- Number of toroidal cuts used to define each coil. Default 5.
%
% Output:
%	rwm.coil(npts,1:3) -- Coil set coordinates in meters. The second index corresponds
%	                      to the x,y, or z coordinate in meters. 
%   rwm.current(npts)  -- Current in each filament. The filaments connecting the RWM 
%					      coils have 0 current, and the last value is irrelevant 
%					     (number of filaments = number of points - 1).
%
%
% Notes:
%	1) Should return geometry aligned to 'NSTX physics' toroidal angle.
%	   0 degrees is at the center of RWM coil 1, and phi increases CCW when
%      viewing NSTX from above, co-Ip, and opposite the direction of increasing
%      RWM coil number.
%	2) Viewed from outside the machine, each 'real' RWM coil current flows CCW for 
%	   a positive value of MDSplus signal "irwm#", thus Br is positive for 
%      positive "irwm#" signal. My coils have the opposite sense, so the taper values
%      are multiplied by -1 below to give the proper field orientation.
% 
% JDL

if length(taper) == 1 
    taper = taper*[1,-1,1,-1,1,-1];
end
if nargin < 2 
    ntorpts = 5;
end
npts = 2*ntorpts + 1;

% Coils are approximately rectangular.
phicens = [0,300,240,180,120,60]*pi/180;  % Center of each coil
phiext = 56*pi/180;						  % Toroidal extent of each coil
R = [1.76,1.76];                          % Major radius of coil
Z = [0.4826,-0.4826];                     % Vertical extent of coil.
nturn = 2;                                % number of electrical turns

coil = zeros(6*npts,3);
current = zeros(6*npts,1);

taper = -nturn*taper;  % Minus sign accounts for the handedness of the coils as I have 
                       % defined them relative to the 'real' orientation described above.
					   % I.e., the returned coils are CW as viewed from outside NSTX.

for i = 0:5
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

