clearvars;

% -------------------------------------------------
% CONSTANTS
% -------------------------------------------------
eps0 = 8.854187817e-12;  % F/m
e = 1.602e-19;
me = 9.11e-31;
mp = 1.672621637e-27;     % kg

% -------------------------------------------------
% BULK PLASMA -------------------------------------
% -------------------------------------------------
Rm = 13;
mi = mp;
ne_m3_bulk = 1e20;
ni_m3_bulk = ne_m3_bulk;
Ti_eV_bulk = 3;
Te_eV_bulk = 10;   % eV

% -------------------------------------------------
% TEST PARTICLES ----------------------------------
% -------------------------------------------------
% T_eV_test = 10;  species_tag = 'e'; m_test = me;
% T_eV_test = 3;  species_tag = 'i'; m_test = mi;
T_eV_test = linspace(0.5,100,100); species_tag = 'e'; m_test = me;
% 
% T_eV_test = linspace(0.5,100,100); species_tag = 'i'; m_test = mi;



% PLASMA


% End plasma


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
theta_deg = asin(1/sqrt(Rm))*180/pi; theta_rad = theta_deg*pi/180;
delta_theta_deg = 90 - asin(1/sqrt(Rm))*180/pi; delta_theta_rad = delta_theta_deg*pi/180;

debye_len = sqrt((Te_eV_bulk./1e3)/(ne_m3_bulk/1e20))*2.35e-5;
lnlam_C = log(9*4/3*pi*ne_m3_bulk*(debye_len).^3);
lnlam_HH = 23.4-1.15*log10(ne_m3_bulk/1e6)+3.45*log10(Te_eV_bulk);

Vth_e = sqrt(2*e*Te_eV_bulk/me);   %%% THERMAL VELOCITY
Vth_i = sqrt(2*e*Ti_eV_bulk/mi);
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

fprintf('\n\n--------------------------------------------------------------\n')
fprintf(' Mirror ratio    : %8.1f\n', Rm)
fprintf(' Loss cone angle : %8.1f (deg) = %8.1f (radian)\n', theta_deg, theta_rad)
fprintf('--------------------------------------------------------------\n')
fprintf(' Te,i (eV) : %8.1f, %8.1f\n', Te_eV_bulk,Ti_eV_bulk)
fprintf(' ne,i (m-3): %8.1e, %8.1e\n', ne_m3_bulk,ni_m3_bulk)
fprintf(' logLam    : %8.1f, %8.1f\n', lnlam_C,lnlam_HH)
fprintf(' Vth_e     : %8.1e\n', Vth_e)
fprintf(' Vth_i     : %8.1e\n', Vth_i)

% Maxwell integrals
nx = 1000;
max_x_eval = 1e5;
min_x_eval = 1e-6;
x = logspace(log10(min_x_eval),log10(max_x_eval),nx);
nt = 1000;
psi = zeros(1,nx);
for ix = 1:nx
    t = logspace(-6,log10(x(ix)),nt);
    psi(ix) = 2/sqrt(pi)*trapz(t,sqrt(t).*exp(-t));
end
psi_prime = deriv(x,psi);
psi_prime2 = 2/sqrt(pi)*sqrt(x).*exp(-x);
fac = 2*(psi + psi_prime - psi/2./x);

figure; hold on; box on;
plot(x,psi,'linewidth',2)
plot(x,psi_prime,'linewidth',2)
plot(x,fac,'linewidth',2)
xlabel('x','fontsize',14)
set(gca,'fontsize',14)
legend('\psi','\psi''','\psi + \psi''-\psi/2x')

Lpar = 1;
ln_lambda = lnlam_C;
for ii = 1:length(T_eV_test)

    V_test = sqrt(2*e*T_eV_test(ii)/m_test);

    x_test_e = V_test^2/Vth_e^2;
    x_test_i = V_test^2/Vth_i^2;

    if x_test_e < min_x_eval
        warning('here')
        psi_xe  = 4*x_test_e^1.5/3/sqrt(pi);
        psip_xe = 2*x_test_e^0.5/sqrt(pi);
    else
        psi_xe = interp1(x,psi,x_test_e);
        psip_xe = interp1(x,psi_prime,x_test_e);
    end

    if x_test_i > max_x_eval
        warning('here')
        psi_xi = 1;
        psip_xi = 0;
    else
        psi_xi = interp1(x,psi,x_test_i);
        psip_xi = interp1(x,psi_prime,x_test_i);
    end

    nu_0_test_e = 4*pi*ne_m3_bulk*e^4*ln_lambda/((4*pi*eps0)^2*m_test^2*V_test^3);

    fac_e = (psi_xe + psip_xe - psi_xe/2/x_test_e);
    fac_i = (psi_xi + psip_xi - psi_xi/2/x_test_i);

    nu_perp_test_e = 2*fac_e*nu_0_test_e;
    nu_perp_test_i = 2*fac_i*nu_0_test_e;
    
    nu_perp = nu_perp_test_e + nu_perp_test_i;
    delta_t(ii) = (delta_theta_deg*pi/180)^2/nu_perp;
    delta_t_par(ii) = Lpar/V_test;
end

figure; hold on; box on;
plot(T_eV_test,delta_t_par,'linewidth',2)
plot(T_eV_test,delta_t,'linewidth',2)
set(gca,'fontsize',14)
legend(strcat('\Deltat_|_|^',species_tag),strcat('\Deltat_{detrap}^',species_tag))
xlabel('T_{test} (eV)','fontsize',14)
set(gca,'yscale','log')

figure; hold on; box on;
plot(T_eV_test,delta_t_par./delta_t,'linewidth',2)
xlabel('T (eV)','fontsize',14)
ylabel('\Deltat_|_|/\Deltat_{detrap}','fontsize',14)
set(gca,'fontsize',14)
% plot(T_test,delta_t_par)



