function [x2,y2,z2,t]=plane_pt_near_pt(a,b,c,d,x,y,z)
% Finds the nearest point on the plane a-d to the point x,y,z
% Based on algorithm by John Burkardt
% JDL 2/24/2012

t = -( a*x + b*y + c*z + d )/( a*a + b*b + c*c );
x2 = x + a*t;
y2 = y + b*t;
z2 = z + c*t;

end