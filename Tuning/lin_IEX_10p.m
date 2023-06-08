function [re, ri, dt, M, t] = lin_IEX_10p(reh, rih)    
% Network with N total inhibition-stabilized excitatory-inhibitory units. 
% M of those pairs we attempt to make "active", and N-M pairs we try to 
% suppress.

%Originally Written by Paul Miller
%Modified and Commented by Connor Zawacki and Yilin Xu

% function called in main, input is initial firing rate for e and i unit
%No WEIX cross connection
%WIEX cross connection
%linear inhibitory unit Firing rate
%linear excitatory unit Firing rate
    
    %clear
    N = 10;             %Total excitatory-inhibitory firing rate unit pairs
    M = 8;              %Total number of active excitatory-inhibitory firing rate unit pairs
    
    dt = 0.001;         %Time step for simulation
    tmax = 100;         %Duration of simulation
    t = 0:dt:tmax;      
    Nt = length(t);
    
    taue = 0.010;       %Time constant for excitatory cells
    taui = 0.010;       %Time constant for inhibitory cells
    
    Wee0 = 3.25;        %Excit.-Excit. self conneciton strength
    Wie0 = -2.4;        %Inhib.-Excit. connection strength
    Wei0 = 5;           %Excit-Inhib. connection strength
    Weix = 0.0;         %Excit-Inhib. cross connection strength
    Wiex = -0.5;        %Inhib.-Excit cross connection strength, Increase cross inhibition will cause increased deviation
    Wii0 = -3;          %Inhib-Inhib self inhibition
    
    I0e = 10;           %Excit. applied current
    I0i = -10;          %Inhib. applied current
    
    %Connection Matricies
    Wee = Wee0*eye(N);
    Wie = Wie0*eye(N);
    Wie = (Wie0-Wiex)*eye(N) + Wiex*ones(N);
    Wei = (Wei0-Weix)*eye(N) + Weix*ones(N);
    Wii = Wii0*eye(N);
    
    Ie = I0e*ones(N,1);
    Ii = I0i*ones(N,1);
    
    %% additional stimulus excitatory current
    stimulus = zeros(N,1);
    stimulus(1:M) = 1;
    Ie_stim = I0e*stimulus;
    Ii_stim = I0i*stimulus;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Inhitial Conditions Calculations
    Wei_tilde = Wei0/(1-Wii0);
    Wie_tilde = (-Wie0)/(1-Wii0);
    Wiex_tilde = (-Wiex)/(1-Wii0);
    I0i_tilde = I0i*Wie_tilde;
    
    jac_stable =  (-Wie0+Wiex)*Wei_tilde - (Wee0-1);
    disp(jac_stable)
    
    ril = I0i/(1-Wii0);
    if ( ril > 0 )
        th_e_tilde = -I0e + (N-M)*Wiex_tilde*I0i;
    else
        th_e_tilde = -I0e;
    end
    
    % Initial firing rate calculated here, default value if no params given
    if nargin<2
        rih = (Wei0*th_e_tilde + I0i*(Wee0-1) )/ ...
            ( (1-Wii0)*(Wee0-1) + Wei0*(Wie0 + (M-1)*Wiex) );
        
        reh = ( th_e_tilde - rih*(Wie0 + (M-1)*Wiex) )/(Wee0-1);
    end

    if ( ril < 0 ) 
        low_rate_input = M*Wiex*rih + I0e
    else
        low_rate_input = M*Wiex*rih + I0e + (Wie0 + (N-M-1)*Wiex)*ril
    end
    
    %Initial firing rates
    re = zeros(N,Nt);
    ri = zeros(N,Nt);
    
    
    %% Alter initial condition here to find more stable states
    %reh = 5;
    %rih = 5;
    %%%%%%%%%
    
    re(1:M,1) = reh;
    ri(1:M,1) = rih;  
    
    sigman = 0.05/sqrt(dt);
    %sigman=0;
    
    for i = 2:Nt
        
        re(:,i) = re(:,i-1) + dt*((Wee-eye(N))*re(:,i-1) + Wie*ri(:,i-1) + Ie + sigman*randn(N,1));
        ri(:,i) = ri(:,i-1) + dt*(Wei*re(:,i-1) + (Wii-eye(N))*ri(:,i-1) + Ii + sigman*randn(N,1));
        
        if (i*dt)>50 && (i*dt)<50.3
            re(:,i) = re(:,i) + Ie_stim*dt;
            ri(:,i) = ri(:,i) + Ii_stim*dt;
        end
        re(:,i) = max(re(:,i),0);
        ri(:,i) = max(ri(:,i),0);
    end
    
    %% processing data here
    %[average1, deviation1] = stable_state_counter(re, dt, M);
    %[average2, deviation2] = stable_state_counter(ri, dt, M);
    
    %disp(deviation1)
    
    
    
end

