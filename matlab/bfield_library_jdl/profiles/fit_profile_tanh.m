function c = fit_profile_tanh(x,y,ystd,quiet,nc,fit_type)
if nargin < 4
    quiet = 0;
end
if nargin < 5
    nc = 9;
end
if nargin < 6
    fit_type = 0;  % 0 is R-Rsep, 1 is psiN
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
if fit_type == 0
    c0(1:4) = [0,0.01,max(y)-min(y),min(y)];
elseif fit_type == 1
    c0(1:4) = [1,0.02,max(y)-min(y),min(y)];
else
    error('Unknown fit_type')
end

% ;1       c0 = SYMMETRY POINT
% ;2       c1 = FULL WIDTH
% ;3       c2 = HEIGHT
% ;4       c3 = OFFSET
% ;5       c4 = SLOPE OR QUADRATIC (IF ZERO DER) INNER
% ;6       c5 = QUADRADIC OR CUBIC (IF ZERO DER) INNER
% ;7       c6 = CUBIC OR QUARTIC (IF ZERO DER) INNER
% ;8       c7 = SLOPE OUTER
% ;9       c8 = QUADRATIC OUTER

clb = []; cub = [];
c=lsqnonlin(@evaluate_tanh_fit_minfun,c0,clb,cub,opts);



    function chi = evaluate_tanh_fit_minfun(c)
        c([6,7,9]) = 0;
        chi = (evaluate_tanh_fit(c,x) - y)./ystd;
    end
end
