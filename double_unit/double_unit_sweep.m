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
N = 2;              %Total excitatory-inhibitory firing rate unit pairs
M = 1;              %Total number of active excitatory-inhibitory firing rate unit pairs           

dt = 0.0001;        %Time step for simulation
tmax = 10;          %Duration of simulation
t = 0:dt:tmax;      %Simulation time vector
Nt = length(t);

taue = 0.010;       %Time constant for excitatory cells
taui = 0.010;       %Time constant for Inhibitory cells

%Wee0 = 1.5;         %Excit.-Excit. connection strength
%Wie0 = -0.4;        %Inhib.-Excit. connection strength
Wei0 = 2.5;         %Excit-Excit connection strength
Weix = 0.5;           %Excit.-Inhibitory Cross connection strength
Wiex = 0;
Wii0 = -1;          %Inhib-Inhib connection strength

I0e = 1;            %Excit. applied current
I0i = 17;           %Inhib. applied current

alpha_e = 1;        %Gain of excit. cells
alpha_i = 0.02;     %Gain of inhib. cells

%Current Matricies
Ie = I0e*ones(N,1);
Ii = I0i*ones(N,1);

sigman = 0.01/sqrt(dt);

%% Set up params for multiple trials
Wee0_vec = 0:0.1:10;
Wie0_vec = 0:0.1:5;
Nvec1 = length(Wee0_vec);
Nvec2 = length(Wie0_vec);
outputmat = zeros(Nvec1, Nvec2);
highlight = [1.6, -0.4];
output_re = zeros(Nvec1*N, Nvec2*Nt);
output_ri = zeros(size(output_re));

%% Simulation
for i = 1:Nvec1
    i
    Wee0 = Wee0_vec(i);
    for j = 1:Nvec2
        Wie0 = -Wie0_vec(j);

        [re, ri] = ...
            linE_QI(N, M, Nt, dt, Wee0, Wie0, Wiex, Wii0, Wei0, Weix, sigman, Ie, Ii, taue, taui, alpha_e, alpha_i);
        outputmat(i,j) = is_bistable(N, M, re(:,ceil(1/dt:end)));

        %output_re = export_fr(N, Nt, re, output_re, i, j);
        %output_ri = export_fr(N, Nt, ri, output_ri, i, j);

        %% Analyze and Plot Highlight Point
        if Wee0 == highlight(1) && Wie0 == highlight(2)
            Wie_tilde = (-Wie0)/(1-Wii0);
            I0i_tilde = I0i*Wie_tilde;
            off_stable = Wee0 - 1 - Wie_tilde*(Wei0-Weix);
            disp("off_stable =" + num2str(off_stable))

            figure(1)
            clf
            sgtitle("Wee0="+num2str(Wee0)+"; Wie0="+num2str(Wie0))
            subplot(2,1,1), 
            plot(t,re), legend("Active","Non-active"), xlabel("Excit."), ylabel("Firing Rate"), 
            subplot(2,1,2)
            plot(t,ri), legend("Active","Non-active"), xlabel("Inhib."), ylabel("Firing Rate")
        end
    end
end

%% Export firing rate matrix re and ri
%writematrix(output_re,'output_re(20,-20).csv')

%writematrix(output_ri, 'output_ri(20,-20).csv')

x = [0, -5];
y = [0, 10];
figure(99), imagesc(x,y,outputmat);
set(gca,'YDir','normal'), xlabel("Wie0"), ylabel("Wee0");

