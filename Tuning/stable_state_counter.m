function [average, deviation] = stable_state_counter(r, dt, M)
    average = zeros(M,1);
    deviation = size(average);

    init_time = round(10/dt);
    init_time = 1;
    for i = 1:M
        average(i) = mean(r(i,init_time:end));
        deviation(i) = std(r(i, init_time:end));
    end

end