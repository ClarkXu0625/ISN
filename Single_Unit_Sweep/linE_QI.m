function [frmat_i,frmat_e] = linE_QI(N, tvec, dt, WEI, WEE, WII, WIE, WEIX, Iapp_i, Iapp_e, theta_i, theta_e, tao_i, tao_e, alpha_e, alpha_i, rmax) 

    frmat_e = zeros(N,numel(tvec));
    Imat_e = frmat_e;
    frmat_i = zeros(N,numel(tvec));
    Imat_i = frmat_i;
    frmat_i(1, 1)= 5;

    %% simulation
    for t = 2:numel(tvec)
        Imat_i(:,t) = WEI*frmat_e(:,t-1) + WII*frmat_i(:,t-1) + Iapp_i(:,t) + WEIX*(sum(frmat_e(:,t-1)) - frmat_e(:,t-1));
        Imat_e(:,t) = WEE*frmat_e(:,t-1) + WIE*frmat_i(:, t-1) + Iapp_e(:,t);
    
        % Linear-E and Quadratic-I
        frmat_e(:,t) = frmat_e(:,t-1) + (dt/tao_e).*(-frmat_e(:,t-1) + alpha_e*(Imat_e(:,t)-theta_e).^2 .*sign(Imat_e(:,t)-theta_e));
        frmat_i(:,t) = frmat_i(:,t-1) + (dt/tao_i).*(-frmat_i(:,t-1) + alpha_i*(Imat_i(:,t)-theta_i).^2 .*sign(Imat_i(:,t)-theta_i));
        
        frmat_i(:,t)=min(frmat_i(:,t),rmax);
        frmat_e(:,t)=min(frmat_e(:,t),rmax);
        frmat_i(:,t)=max(frmat_i(:,t),0);
        frmat_e(:,t)=max(frmat_e(:,t),0);       
    end
end