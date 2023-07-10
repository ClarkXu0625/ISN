% Given the number of unit, firing rate of excit. and inhib. units, return
% if the units showing multi-stablity, or firing at uniform firing rate
function [bistable] = is_bistable(N, re, ri)
    excit_average = zeros(1,N);
    inhib_average = zeros(1,N);
    for i = 1:N
        excit_average(i) = mean(re(i,:));
        inhib_average(i) = mean(ri(i,:));
    end
    disp(std(excit_average))
    disp(std(inhib_average))
    bistable = (std(excit_average)>1 || std(inhib_average)>1);
end