% Network with N total inhibition-stabilized excitatory-inhibitory units. 
% M of those pairs we attempt to make "active", and N-M pairs we try to 
% suppress.

%In particular, this simulation displays bistability within a single 
%unit given a quadratic firing rate equation for the inhibitory units

%No WEIX cross connection
%No WIEX cross connection
%nonlinear inhibitory unit Firing rate
%linear excitatory unit Firing rate

% Requirement: current working directory should be double_unit

clear

%% Add all subfolders of current directory into matlab session search
current_path = pwd;
addpath(genpath(current_path))

%% Parameters
N = 5;              %Total excitatory-inhibitory firing rate unit pairs

%M = 2;              %Total number of active excitatory-inhibitory firing rate unit pairs           

dt = 0.0001;        %Time step for simulation
tmax = 10;          %Duration of simulation
t = 0:dt:tmax;      %Simulation time vector
Nt = length(t);

taue = 0.010;       %Time constant for excitatory cells
taui = 0.010;       %Time constant for Inhibitory cells

%Wee0 = 1.5;         %Excit.-Excit. connection strength
%Wie0 = -0.4;        %Inhib.-Excit. connection strength
Wei0 = 2.5;          %Excit-Excit connection strength
Weix = 0.5;          %Excit.-Inhibitory Cross connection strength
Wiex = 0;
Wii0 = -1;          %Inhib-Inhib connection strength

I0e = 1;            %Excit. applied current
I0i = 17;           %Inhib. applied current

alpha_e = 1;        %Gain of excit. cells
alpha_i = 0.02;     %Gain of inhib. cells

%Current Matricies
Ie = I0e*ones(N,1);
Ii = I0i*ones(N,1);

sigman = 0.005/sqrt(dt);

%% Set up params for multiple trials
max_vec1 = 6;   % maximum value for vector 1
max_vec2 = 5;   % maximum value for vector 2
Wee0_vec = 0:0.2:max_vec1;
Wie0_vec = 0:0.2:max_vec2;
Nvec1 = length(Wee0_vec);
Nvec2 = length(Wie0_vec);
outputmat = zeros(N+1, Nvec1, Nvec2);
stable_state = zeros(Nvec1, Nvec1);
highlight = [-1, -1];

for M = 0:N    
    %% Simulation
    for i = 1:Nvec1
        progress = M*Nvec1+i;
        progress

        Wee0 = Wee0_vec(i);

        for j = 1:Nvec2
            Wie0 = -Wie0_vec(j);
    
            [re, ri] = ...
                linE_QI(N, M, Nt, dt, Wee0, Wie0, Wiex, Wii0, Wei0, Weix, sigman, Ie, Ii, taue, taui, alpha_e, alpha_i);
            [bistable, Nstate] = is_bistable(N, M, re(:,ceil(1.5/dt:end)));
            outputmat(M+1,i,j) = bistable;
            stable_state(i,j) = stable_state(i,j) + Nstate;
        end
    end
end
%% Export stable state to a csv file
output = export_bistable(N, Nvec1, Nvec2, outputmat);
writematrix(output,'5ss(6,-5).csv')

x = [0, -max_vec2];
y = [0, max_vec1];
figure(99), imagesc(x,y,stable_state);
set(gca,'YDir','normal'), xlabel("Wie0"), ylabel("Wee0");
title("Num of stabel states, Wie0-Wee0 sweep, Wei0 = " + num2str(Wei0))
c = colorbar;
c.Label.String = "stable state number";
