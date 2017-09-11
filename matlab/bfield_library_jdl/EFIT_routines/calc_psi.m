function [psi] = calc_psi(g,R1,Z1,quiet)
%[psi] = calc_psi(g,R1,Z1,quiet)
% R1, Z1 can be vectors
% quiet optional  
%--- Just applies g.ip_sign to get_psi_bicub!
if nargin < 4
    quiet = 0;
end
psi = g.ip_sign*get_psi_bicub(g,R1,Z1,quiet);