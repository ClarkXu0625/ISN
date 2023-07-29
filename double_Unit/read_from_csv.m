N = 5;
%% Add all subfolders of current directory into matlab session search
current_path = pwd;
addpath(genpath(current_path))

% read matrix from csv
[output, Nvec1, Nvec2] = read_bistable("5ss(6,-5).csv", N);
stable_states = zeros(Nvec1, Nvec2);

for M = 0:N
    Nstate = factorial(N)/(factorial(N-M)*factorial(M));
    disp(Nstate)
    stable_states = stable_states + Nstate*reshape(output(M+1,:,:), [Nvec1, Nvec2]);
end

x = [0, -5];
y = [0, 6];
figure(99), imagesc(x,y,stable_states);
set(gca,'YDir','normal'), xlabel("Wie0"), ylabel("Wee0");
title("Num of stabel states, Wie0-Wee0 sweep, Wei0 = " + num2str(Wei0)+"; N = "+num2str(N))
c = colorbar;
c.Label.String = "stable state number";

% Function read the csv from given file name. The csv file contains the
% bistable information for a system of N unit.
function [output, Nvec1, Nvec2] = read_bistable(filename, N)
    
    A = readmatrix(filename);
    Nvec1 = size(A, 1)/(N+1);
    Nvec2 = size(A, 2);
    output = zeros(N, Nvec1, Nvec2);
    for i = 1: N+1
        output(i,:,:) = A((i-1)*Nvec1+1:i*Nvec1, :);
    end
end
