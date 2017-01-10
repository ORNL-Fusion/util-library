clearvars;
% !  pcurr_type       character, specifies the parameterization of the profile
% !                     | - X for parametrization of I-prime(s), _ for I(s)
% !    gauss_trunc      X    Truncated Gaussian - I-prime
% !    two_power        X    Two powers - ac(0) * (1 - s ** ac(1)) ** ac(2)
% !    two_power_gs     X    Two powers with gaussian peaks -
% !                          ac(0) * ((1 - s ** ac(1)) ** ac(2))*(1 + Sum[ac(i)*Exp(-(s - ac(i+1))/ac(i+2)) ** 2])
% !    sum_atan         _    sum of arctangents
% !    power_series_I   _    Power series for I(s) (NOT default)
% !    Akima_spline_Ip  X    Akima spline for I-prime(s)
% !    Akima_spline_I   _    Akima spline for I(s)
% !    cubic_spline_Ip  X    cubic spline for I-prime(s)
% !    cubic_spline_I   _    cubic spline for I(s)
% !    pedestal         _    Pedestal profile
% !    rational         _    Rational function (ratio of polynomials)
% !    line_segment_Ip  X    Line segments for I-prime(s)
% !    line_segment_I   _    Line segments for I(s)
% !    power_series     X    Power Series for I-prime(s) (Default)


s = 0:0.01:1;

ac = [.85,2,2];

% Te_core = 3000;
% Ti_core = 1000;
% n_core = 1.e20;

Te_core = 5000;
Ti_core = 3500;
n_core = 1.e20;

e = 1.6e-19;

p_core = n_core*e*(Te_core + Ti_core);

prof = ac(1)*((1 - s.^ac(2)).^ac(3));
prof2 = (1-s.^4).*(1-s);

figure; hold on; box on;
plot(s,p_core*prof,'b','linewidth',2)
plot(s,p_core*prof2,'k','linewidth',1)
xlabel('s','fontsize',14,'fontweight','bold')
ylabel('P [Pa]','fontsize',14,'fontweight','bold')
set(gca,'fontsize',14,'fontweight','bold')

B0 = 2.5;
mu0 = 4*pi*1e-7;
beta = p_core*prof*2*mu0/B0.^2;
beta2 = p_core*prof2*2*mu0/B0.^2;


figure; hold on; box on;
plot(s,beta*100,'b','linewidth',2)
plot(s,beta2*100,'k','linewidth',1)
xlabel('s','fontsize',14,'fontweight','bold')
ylabel('\beta [%]','fontsize',14,'fontweight','bold')
set(gca,'fontsize',14,'fontweight','bold')