

[re, ri, dt, M, t] = lin_IEX_10p();
[average1, deviation1] = stable_state_counter(re, dt, M);
    %[average2, deviation2] = stable_state_counter(ri, dt, M);
    
    %disp(deviation1)

%% plot data
figure(1)
clf
subplot(2,1,1)
plot(t,re)
title("E-unit")
xlabel("time (s)")
ylabel("firing rate")

subplot(2,1,2)
plot(t,ri)
title("I-unit")
xlabel("time (s)")
ylabel("firing rate")