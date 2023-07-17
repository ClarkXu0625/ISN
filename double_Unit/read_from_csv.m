N = 2;
[output, Nvec1, Nvec2] = read_bistable("stablestates(6,-5).csv", N);
stable_states = zeros(Nvec1, Nvec2);

for M = 0:2
    Nstate = factorial(N)/(factorial(N-M)*factorial(M))*bistable;
    disp(Nstate)
    stable_states = stable_states + Nstate*reshape(output(M+1,:,:), [Nvec1, Nvec2]);
end

figure(99), imagesc(stable_states);
set(gca,'YDir','normal'), xlabel("Wie0"), ylabel("Wee0");
title("Num of stabel states, Wie0-Wee0 sweep, Wei0 = " + num2str(Wei0))
c = colorbar;
c.Label.String = "stable state number";


function [output, Nvec1, Nvec2] = read_bistable(filename, N)
    
    A = readmatrix(filename);
    Nvec1 = size(A, 1)/(N+1);
    Nvec2 = size(A, 2);
    output = zeros(N, Nvec1, Nvec2);
    for i = 1: N+1
        output(i,:,:) = A((i-1)*Nvec1+1:i*Nvec1, :);
    end
end
