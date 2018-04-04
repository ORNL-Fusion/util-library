function [Rs,Zs,Ps] = symmetrize_and_move_to_hfp(R,Z,P,nfp,stellsym,fp_want,move_to_2nd_hfp)
    if fp_want <= 0 || fp_want > nfp
        error('bad fp want: %d',fp_want)
    end
    if move_to_2nd_hfp && ~stellsym
        error('Must have stellsym if move_to_2nd_hfp')
    end
    Rs = zeros(size(R));
    Zs = Rs;
    Ps = Rs;
    for i = 1:numel(R)
        RZP = [R(i),Z(i),P(i)];
        rpz = symmetrize([RZP(1),RZP(3),RZP(2)],nfp,stellsym);
        Rs(i) = rpz(1);
        Ps(i) = rpz(2);
        Zs(i) = rpz(3);
    end
    if move_to_2nd_hfp
        Ps = -Ps + 2*pi/nfp;
        Zs = -Zs;
    end
    
    % Now everything is in first FP (or HFP if stellsym)
    Ps = Ps + (fp_want-1)/nfp*2*pi;
    
end
