clearvars;

data = return_AMJUEL_data;

mp = 1.67e-27;
matom = 2*mp;  % 2, 184
mmolc = 2*matom;
Tatom = 3;   % 3
Tatom_slow = 0.026;
Tmolc = 0.026;

vatom = sqrt(2*Tatom*1.6e-19/matom);
vatom_slow = sqrt(2*Tatom_slow*1.6e-19/matom);
vmolc = sqrt(2*Tmolc*1.6e-19/mmolc);

% rH0 = 5.3e-11;
rH2 = 1.37e-10;
n0_coll_xs = pi*(2*rH2)^2;
nn0 = 3.3e19;

% suffix = '50'; element='w';
% suffix = '96'; element='h';
suffix = '12'; element='h';
fname = ['C:\Work\ADAS\adf11_all\scd',suffix,'\','scd',suffix,'_',element,'.dat'];  % Effective ionization coefficients (cm^-3/s)
scd = read_adas_adf11_file(fname); scd1 = squeeze(scd.coeff(:,:,1));
fname = ['C:\Work\ADAS\adf11_all\acd',suffix,'\','acd',suffix,'_',element,'.dat'];  % Effective recombination coefficients (cm^-3/s)
acd = read_adas_adf11_file(fname); acd1 = squeeze(acd.coeff(:,:,1));

if 1
    add_slow = 1;

    n0_mfp = 1/(sqrt(2)*n0_coll_xs*nn0)
    
    % Te_test = logspace(log10(1),4,20);
    % ne_test = logspace(8,16,5);
    %
    Te_test = linspace(.1,20,100);
    ne_test = logspace(log10(1e13),log10(1e13),1);  % cm3
    
    for i=1:length(ne_test)
        sv_iz_adas(:,i) = interp_adas_rate_coefficient(Te_test,ne_test,scd.te,scd.ne,scd1);
        sv_iz(:,i) = eval_AMJUEL_H4_fit(data.H4.reaction_215,ne_test(i),Te_test);               
        sv_rc(:,i) = eval_AMJUEL_H4_fit(data.H4.reaction_218,ne_test(i),Te_test);
        sv_rc_adas(:,i) = interp_adas_rate_coefficient(Te_test,ne_test,acd.te,acd.ne,acd1);
        sv_H2_diss_iz(:,i) = ...
            eval_AMJUEL_H4_fit(data.H4.reaction_225 ,ne_test(i),Te_test) + ...
            eval_AMJUEL_H4_fit(data.H4.reaction_229 ,ne_test(i),Te_test) + ...
            eval_AMJUEL_H4_fit(data.H4.reaction_2210,ne_test(i),Te_test);
    end
    sv_cx = eval_AMJUEL_H3_fit(data.H3.reaction_318,Tatom,Te_test);
    
    
    s = styflipper(length(ne_test));
    
    if 0
        % Rates for low Te
        figure; hold on; box on
        for i=1:length(ne_test)
            plot(Te_test,sv_iz(:,i),'r','linewidth',2,'linestyle',char(s{i}))
            plot(Te_test,sv_rc(:,i),'b','linewidth',2,'linestyle',char(s{i}))
            plot(Te_test,sv_H2_diss_iz(:,i),'k','linewidth',2,'linestyle',char(s{i}))
        end
        plot(Te_test,sv_cx,'m','linewidth',2,'linestyle',char(s{i}))
        % set(gca,'xscale','log')
        set(gca,'yscale','log')
        xlabel('T_e_,_i (eV)','fontsize',12)
        ylabel('<\sigmav> (cm^3/s)','fontsize',12)
        set(gca,'fontsize',12)
        axis([.1,10,10^-16,10^-7])
        legend('D iz','D rc (rad+3bdy)','D2 iz+ds','D cx')
    end
    
    
    for i=1:length(ne_test)
        mfp_iz_atom(:,i) = vatom./(ne_test(i)*sv_iz(:,i));
        mfp_iz_atom_adas(:,i) = vatom./(ne_test(i)*sv_iz_adas(:,i));
        mfp_iz_atom_slow(:,i) = vatom_slow./(ne_test(i)*sv_iz(:,i));
        mfp_rc_atom(:,i) = vatom./(ne_test(i)*sv_rc(:,i));
        mfp_rc_atom_adas(:,i) = vatom./(ne_test(i)*sv_rc_adas(:,i));
        mfp_rc_atom_slow(:,i) = vatom_slow./(ne_test(i)*sv_rc(:,i));
        mfp_molc(:,i) = vmolc./(ne_test(i)*sv_H2_diss_iz(:,i));
        mfp_cx(:,i) = vatom./(ne_test(i)*sv_cx);
    end
    
    
    figure; hold on; box on
    for i=1:length(ne_test)
        plot(Te_test,mfp_iz_atom(:,i),'r','linewidth',2,'linestyle',char(s{i}))
%         plot(Te_test,mfp_iz_atom_adas(:,i),'c','linewidth',2,'linestyle',char(s{i}))
        plot(Te_test,mfp_rc_atom(:,i),'b','linewidth',2,'linestyle',char(s{i}))
%         plot(Te_test,mfp_rc_atom_adas(:,i),'y','linewidth',2,'linestyle',char(s{i}))
        plot(Te_test,mfp_molc(:,i),'k','linewidth',2,'linestyle',char(s{i}))
        if add_slow
            plot(Te_test,mfp_iz_atom_slow(:,i),'r','linewidth',2,'linestyle','--')
            plot(Te_test,mfp_rc_atom_slow(:,i),'b','linewidth',2,'linestyle','--')
        end
        plot(Te_test,mfp_cx(:,i),'m','linewidth',2,'linestyle',char(s{i}))
    end
    
    plot(Te_test,n0_mfp*ones(size(Te_test)),'k--','linewidth',2)
    % set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('T_e_,_i (eV)','fontsize',12)
    ylabel('\lambda_{mfp} (m)','fontsize',12)
    set(gca,'fontsize',12)
    axis([.1,20,10^-3,10^3])
    if add_slow
        h = legend('D iz','D rc (rad+3bdy)','D2 iz+ds','D iz (slow)','D rc slow (rad+3bdy)','D cx','D2-D2');
    else
        h = legend('D iz','D rc (rad+3bdy)','D2 iz+ds','D cx','D2-D2');
    end
    set(h,'fontsize',12)
end




% Rates for two temperature ranges
if 0
    clear Te_test ne_test;
    Te_test = logspace(-1,3,1000);
    ne_test = logspace(12,14,3);  % cm3
    s = styflipper(length(ne_test));
    
    for i=1:length(ne_test)
        sv_iz(:,i) = eval_AMJUEL_H4_fit(data.H4.reaction_215,ne_test(i),Te_test);
        sv_rc(:,i) = eval_AMJUEL_H4_fit(data.H4.reaction_218,ne_test(i),Te_test);
        sv_H2_diss_iz(:,i) = ...
            eval_AMJUEL_H4_fit(data.H4.reaction_225 ,ne_test(i),Te_test) + ...
            eval_AMJUEL_H4_fit(data.H4.reaction_229 ,ne_test(i),Te_test) + ...
            eval_AMJUEL_H4_fit(data.H4.reaction_2210,ne_test(i),Te_test);
    end
    sv_cx = eval_AMJUEL_H3_fit(data.H3.reaction_318,Tatom,Te_test);
    figure; hold on; box on
    for i=1:length(ne_test)
        plot(Te_test,sv_iz(:,i),'r','linewidth',2,'linestyle',char(s{i}))
        plot(Te_test,sv_rc(:,i),'b','linewidth',2,'linestyle',char(s{i}))
        plot(Te_test,sv_H2_diss_iz(:,i),'k','linewidth',2,'linestyle',char(s{i}))
    end
    plot(Te_test,sv_cx,'m','linewidth',2,'linestyle',char(s{i}))
    % set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('T_e_,_i (eV)','fontsize',12)
    ylabel('<\sigmav> (cm^3/s)','fontsize',12)
    set(gca,'fontsize',12)
    axis([.1,10,10^-16,10^-7])
    legend('D iz','D rc (rad+3bdy)','D2 iz+ds','D cx')
    
    
    figure; hold on; box on
    for i=1:length(ne_test)
        plot(Te_test,sv_iz(:,i),'r','linewidth',2,'linestyle',char(s{i}))
        plot(Te_test,sv_rc(:,i),'b','linewidth',2,'linestyle',char(s{i}))
        plot(Te_test,sv_H2_diss_iz(:,i),'k','linewidth',2,'linestyle',char(s{i}))
    end
    plot(Te_test,sv_cx,'m','linewidth',2,'linestyle',char(s{i}))
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('T_e_,_i (eV)','fontsize',12)
    ylabel('<\sigmav> (cm^3/s)','fontsize',12)
    set(gca,'fontsize',12)
    axis([1,1000,10^-16,10^-7])
    legend('D iz','D rc (rad+3bdy)','D2 iz+ds','D cx')
    
end
