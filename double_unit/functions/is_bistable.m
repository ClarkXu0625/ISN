% Given the number of unit, firing rate of excit. and inhib. units, return
% if the units showing multi-stablity, or firing at uniform firing rate
% Parameter N is number of units, and M is number of "on" unit
function [bistable] = is_bistable(N, M, re)
    bistable = 1;

    % Ensure firing rate reaches a stable state, not oscillating
    for i = 1:M
        if std(re(i,:))>0.5
            bistable = 0;
        end
    end

    if M == N   % All units are on
        for i = 1:M
            if mean(re(i,:)) < 0.5
                bistable = 0;
            end
            
        end
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
            if mean(re(i,:)) < 1
                bistable = 0;
            end
        end

    end
end