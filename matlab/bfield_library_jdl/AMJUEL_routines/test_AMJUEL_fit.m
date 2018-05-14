clearvars;

check_Ralph_plot = 0;
check_cx_data = 0; 
sanity_check_H3_323 = 0;  % Try to reproduce AMJUEL plot

data = return_AMJUEL_data;

check_236A0 = 1;
if check_236A0
    Te_test = logspace(0,4,25);
    ne_test = [1e10,1e12,1e13,1e14,1e15];
    
    for i=1:length(ne_test)
        sv_rc(:,i) = eval_AMJUEL_H4_fit(data.H4.reaction_236A0,ne_test(i),Te_test);
    end
    
    s = styflipper(length(ne_test));
    syms = ['o','x','.','>','s'];
    figure; hold on; box on
    for i=1:length(ne_test)
        plot(Te_test,sv_rc(:,i),'b','linewidth',2,'linestyle',char(s{i}),'marker',syms(i))
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('T_e (eV)','fontsize',12)
    ylabel('<\sigmav> (cm^3/s)','fontsize',12)
    set(gca,'fontsize',12)
    axis([1e0,1e4,10^-16,10^-11])
    grid on;
    
    suffix = '93'; element='c';
    fname = ['C:\Work\ADAS\adf11_all\acd',suffix,'\','acd',suffix,'_',element,'.dat'];  % Effective recombination coefficients (cm^-3/s)
    acd = read_adas_adf11_file(fname);
    charge_index = 1;
    for i = 1:length(ne_test)
        sv_rc_adas(:,i) = interp_adas_rate_coefficient(Te_test,ne_test(i),acd.te,acd.ne,acd.coeff(:,:,charge_index));
    end
    for i=1:length(ne_test)
        plot(Te_test,sv_rc_adas(:,i),'r','linewidth',2,'linestyle',char(s{i}),'marker',syms(i))
    end
%     aa=1;
end

if sanity_check_H3_323
    Ti_test = logspace(-1,3,100);
    E0_test = logspace(-1,2,5);
    
    for i=1:length(E0_test)
        sv(:,i) = eval_AMJUEL_H3_fit(data.H3.reaction_323,E0_test(i),Ti_test);
    end
    
    s = styflipper(length(E0_test));
    figure; hold on; box on
    for i=1:length(E0_test)
        plot(Ti_test,sv(:,i),'r','linewidth',2,'linestyle',char(s{i}))
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('T_i (eV)','fontsize',12)
    ylabel('<\sigmav> (cm^3/s)','fontsize',12)
    set(gca,'fontsize',12)
    axis([0.1,1000,10^-14,1e-5])
    
end


if check_cx_data
%     Ti_test = logspace(-1,3,100);
    Ti_test = logspace(0,3,100);
    E0_test = logspace(-1,2,5);
%     E0_test = 0.025;
    
    for i=1:length(E0_test)
        sv_cx(:,i) = eval_AMJUEL_H3_fit(data.H3.reaction_318,E0_test(i),Ti_test);
    end
    
    s = styflipper(length(E0_test));
    figure; hold on; box on
    for i=1:length(E0_test)
        plot(Ti_test,sv_cx(:,i)*1e-6,'r','linewidth',2,'linestyle',char(s{i}))
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('T_i (eV)','fontsize',12)
    ylabel('<\sigmav> (m^3/s)','fontsize',12)
    set(gca,'fontsize',12)
    axis([1,1000,10^-16,2e-13])
    
end



if check_Ralph_plot
    % Te_test = logspace(log10(1),4,20);
    % ne_test = logspace(8,16,5);
    %
    Te_test = linspace(0.1,5,100);
    ne_test = logspace(13,15,3);
    
    for i=1:length(ne_test)
        sv_iz(:,i) = eval_AMJUEL_H4_fit(data.H4.reaction_215,ne_test(i),Te_test);
        sv_rc(:,i) = eval_AMJUEL_H4_fit(data.H4.reaction_218,ne_test(i),Te_test);
    end
    
    
    s = styflipper(length(ne_test));
    figure; hold on; box on
    for i=1:length(ne_test)
        plot(Te_test,sv_iz(:,i),'r','linewidth',2,'linestyle',char(s{i}))
        plot(Te_test,sv_rc(:,i),'b','linewidth',2,'linestyle',char(s{i}))
    end
    % set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('T_e (eV)','fontsize',12)
    ylabel('<\sigmav> (cm^3/s)','fontsize',12)
    set(gca,'fontsize',12)
    % axis([1,2e4,10^-12,10^-7])
    axis([0,5,10^-16,10^-8])
end