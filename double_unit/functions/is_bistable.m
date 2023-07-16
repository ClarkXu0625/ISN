% Given the number of unit, firing rate of excit. and inhib. units, return
% if the units showing multi-stablity, or firing at uniform firing rate
% Parameter N is number of units, and M is number of "on" unit
% Given firing rate vector in shape of (N x Nt)

% Decide whether the given network is showing multistability in given M and
% N. If all units are on, all neuron should have firing rate above zero,
% not in quiescence state. If all neurons are off, they should be in
% quiescence state. Otherwise, column 1:M should be on, and column M+1:N 
% should be off.

% If the given firing rate vector fulfills all the requirements, the given
% network should at least have N!/[(N-M)!*M!] different stable state,
% returned as Nstate.
% To evaluate the multistability in the network for downstream tasks, this
% function would also return bistable (0 or 1) to help analyze
% multistability patterns.

function [bistable, Nstate] = is_bistable(N, M, re)
    bistable = 1;

    if M == N   % All units are on

    elseif M == 0  % All units are off
        
    else
        
        % "off" units are at low firing rate
        for i = M+1:N
            if mean(re(i,:)) > 0.5
                bistable = 0;
            end
        end

        % "on" unit should not stay in quescence status
        for i = 1:M
            if mean(re(i,:)) < 1 || mean(re(i,:)) > 90
                bistable = 0;
            end            
        end
    end
    Nstate = factorial(n)/(factorial(N-M)*factorial(M))*bistable;
end