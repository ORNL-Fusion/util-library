function [A,B,C,D] = plane_exp2imp_3d(p1,p2,p3)
% Converts an explicit plane (3 points) to implicit form in 3D.
% Implicit form of a plane in 3D is:
% A*X + B*Y + C*Z + D = 0
% Based on a routine by John Burkardt
% JDL 2/28/2011

N = cross(p2-p1,p3-p1);
A=N(1);
B=N(2);
C=N(3);
D = -dot(p2,N);