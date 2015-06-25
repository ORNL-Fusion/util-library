clearvars;
% function point_inside_hexahedron(P,H)
% P is point to be evaluated P = [X,Y,Z], in meters;
% H defines the hexahedron. Faces assumed to be planar
% !           8----7  <-- Second toroidal face (it+1)
% !          /|   /|
% !    Z    / 5--/-6
% !    ^   4----3 /
% !    |   |/   |/
% !    |   1----2  <-- first toridal face (it)
%      --> R
% H(#,1:3) = [X,Y,Z] in meters
% Plan: Define each face, define normal, evaluate height with positive
% 'inside', if all heights are positive then point is inside


H(1,:) = [0,0,0];
H(2,:) = [1,0,0];
H(3,:) = [1,0,1];
H(4,:) = [0,0,1];
H(5,:) = [0,1,0];
H(6,:) = [1,1,0];
H(7,:) = [1,1,1];
H(8,:) = [0,1,1];

figure; hold on;
plot3(H(:,1),H(:,2),H(:,3),'ko-')
