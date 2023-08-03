% Given the number of unit, firing rate of excit. and inhib. units, return
% if the units showing multi-stablity, or firing at uniform firing rate
% Parameter N is number of units, and M is number of "on" unit
% Given firing rate vector in shape of (N x Nt), and should be the firing
% rate vector of an attractor network that have already entered the stable
% state.

% First, the given network should reach a stable state, not oscillating.
% Them, Decide whether the given network is showing multistability in given 
% M and N. 
% If all units are on, all neuron should have firing rate above zero,
% not in quiescence state. If all neurons are off, they should be in
% quiescence state. Otherwise, column 1:M should be on, and column M+1:N 
% should be off.

% If the given firing rate vector fulfills all the requirements, the given
% network should at least have N!/[(N-M)!*M!] different stable state,
% returned as Nstate.
% To evaluate the multistability in the network for downstream tasks, this
% function would also return bistable (0 or 1) to help analyze
% multistability patterns.

function [bistable, Nstate] = is_bistable(N, M, re, dt, start_t)
    bistable = 1;

    % Ensure firing rate reaches a stable state, not oscillating
    for i = 1:M
        %disp(std(re(i , floor(3/dt):floor(start_t/dt-1))))
        if std(re(i , floor(3/dt):floor(start_t/dt-1))) > 0.01 || ... % no oscillation
            mean(re(i,:)) > 90  % not firing at rmax
            bistable = 0;

        end
        
        % Unstable fixed point; if different firing rate before and after
        % noise is introduced, not in an stable fixed point
        if abs(re(i,floor(4.9/dt)) - re(i,floor(9/dt))) > 0.5
            bistable = 0;
            %disp("8")
        end
    end

    if M == N   % All units are on, all unit should fire above 0.5Hz
        for i = 1:M
            if mean(re(i,:)) < 0.5
                bistable = 0;
            end            
        end

    elseif M == 0  % All units are off, all units should below 1Hz.
        for i = 1:M
            if mean(re(i,:)) > 1
                bistable = 0;
            end            
        end
    else   % Units are partially turned on       
        % "off" units are at low firing rate

        for i = M+1:N
            if mean(re(i, floor(3/dt):floor(start_t/dt)-1)) > 0.5
                bistable = 0;
            end
        end

        % "on" unit should not stay in quescence status, and should not
        % always fire at rmax
        for i = 1:M
            period = re(i, floor(3/dt):floor(start_t/dt)-1);
            if mean(period) < 1 || mean(period) > 90
                bistable = 0;
            end            
        end
    end
    Nstate = factorial(N)/(factorial(N-M)*factorial(M))*bistable;
end