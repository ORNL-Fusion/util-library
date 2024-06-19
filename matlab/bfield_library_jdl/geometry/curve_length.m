function L = curve_length(r,z,NORM)
% Returns length along curve using linear line segments
% L goes from 0 to the total length
% If NORM then L goes from 0 to 1

if nargin < 3
    NORM = false;
end

L = [0; cumsum(sqrt(diff(r).^2  + diff(z).^2))];

if NORM
    L = L./L(end);
end