function [r,phi,z]=line_follow_xdr(r,phi,z,distance,num_phi,bfield_xdr)

% r = 6.0;
% phi = 0;
% z = 0;
% distance = pi;
% num_phi = 100;
% 
% out_path = 'C:\Work\Stellarator\ALL_W7X_WORK\xdr_dump_read\OUTPUT\';
% fname = 'field181x181x96.w7x.1000_1000_1000_1000_+0750_+0750.vac.out';
% bfield_xdr = read_xdr_dump_file(out_path,fname);

% Follows a field line starting at r,phi,z.  Tries to follow it for a phi
% angle distance (i.e. to do one toroidal transit, set distance to 2pi).
% If the field line hits the wall,  the integration is stopped and the flag
% ithit is set to one.  Outputs the arrays r,phi,z, which are points along 
% the field line
%--------------------------------------------------------------------------


%if the phi points are not specified just specify initial and final points
% of the integration
% if nargin<7
%     phi_points=[phi phi+distance];
% else
    phi_points=linspace(phi,phi+distance,num_phi);
% end


%ODE45 settings
% options = odeset('RelTol',1e-9,'AbsTol',1e-10);
options = odeset('RelTol',1e-3,'AbsTol',1e-6,'initialstep',0.5*pi/180,'maxstep',1.0*pi/180);
X0=[r z]; %initial conditions

[phi,X] = ode45(@field_line_derivs,phi_points,X0,options,bfield_xdr);  

r = X(:,1);
z = X(:,2);


function [dXdphi] = field_line_derivs(phi,X,bfield)
%this is the derivative function for field line integration

%pull components from solution vector
r = X(1);
z = X(2);
% phi/(2*pi)

% P_x = r*cos(phi);
% P_y = r*sin(phi);
% P_z = z;

% Get the field at r,phi,z
% [Bx,By,Bz] = hsxfield_grid_interp(r,phi,z,current,taper);
% [Bx,By,Bz]=calc_b_bs_JL(r,phi,z,current,taper);
% [Bx,By,Bz]=bfield_bs_jdl(P_x,P_y,P_z,coil,current);
% [Bx,By,Bz]=bfield_bs_jc(P_x,P_y,P_z,coil,current);

[bval,idiv] = bint_xdr(bfield,[r,phi,z]);
% [bval,idiv] = bint_xdr_interp(bfield,[r,phi,z]);

if idiv ~= 0
    bval = [0,1,0];
end

% Convert to cyl coords
% Br = Bx*cos(phi) + By*sin(phi);
% Bphi = -Bx*sin(phi) + By*cos(phi);

Br = bval(1);
Bphi = bval(2);
Bz = bval(3);

% Return derivs
drdphi = (r*Br)/Bphi;
dzdphi = (r*Bz)/Bphi;
dXdphi = [drdphi dzdphi]';




