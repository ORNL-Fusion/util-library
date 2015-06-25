function c = fit_profile_tanh(x,y,ystd,quiet,nc)
if nargin < 4
    quiet = 0;
end
if nargin < 5
    nc = 9;
end
tolfun = 1e-12;
tolx = 1e-12;

if quiet
    qval = 'off';
else
    qval = 'final';
end
opts=optimoptions('lsqnonlin','TolFun',tolfun,'TolX',tolx,'Display',qval);

c0 = zeros(1,nc);
c0(1:4) = [0,0.01,max(y)-min(y),min(y)];

% ;       c0 = SYMMETRY POINT
% ;       c1 = FULL WIDTH
% ;       c2 = HEIGHT
% ;       c3 = OFFSET
% ;       c4 = SLOPE OR QUADRATIC (IF ZERO DER) INNER
% ;       c5 = QUADRADIC OR CUBIC (IF ZERO DER) INNER
% ;       c6 = CUBIC OR QUARTIC (IF ZERO DER) INNER
% ;       c7 = SLOPE OUTER
% ;       c8 = QUADRATIC OUTER

clb = []; cub = [];
c=lsqnonlin(@evaluate_tanh_fit_minfun,c0,clb,cub,opts);



    function chi = evaluate_tanh_fit_minfun(c)
        chi = (evaluate_tanh_fit(c,x) - y)./ystd;
    end
end
