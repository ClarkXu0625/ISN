%% Add all subfolders of current directory into matlab session search
current_path = pwd;
addpath(genpath(current_path))

%% Read from csv
output_re = readmatrix('output_re(20,-20).csv');
output_ri = readmatrix('output_ri(20,-20).csv');

N = 2;
Wee0_vec = 0:5:20;
Wie0_vec = 0:5:20;
Nvec1 = length(Wee0_vec);
Nvec2 = length(Wie0_vec);

outputmat = zeros(Nvec1, Nvec2);
for i = 1:Nvec1
    for j = Nvec2
        re_vec = read_fr(N, Nt, output_re, i, j);
        ri_vec = read_fr(N, Nt, output_ri, i, j);
        outputmat(i,j) = is_bistable(N, M, re, ri);
    end
end

x = [0, -20];
y = [0, 20];
figure(99), imagesc(x,y,outputmat);
set(gca,'YDir','normal'), xlabel("WEI"), ylabel("WEE");