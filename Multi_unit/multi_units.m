% Network with N total inhibition-stabilized excitatory-inhibitory units. 
% M of those pairs we attempt to make "active", and N-M pairs we try to 
% suppress.

% In particular, this simulation displays bistability within a single 
% unit given a quadratic firing rate equation for the inhibitory units

% Nonlinear inhibitory unit Firing rate
% Linear excitatory unit Firing rate

% Noise duration could be manually set

clear
N = 2;              %Total excitatory-inhibitory firing rate unit pairs
M = 1;              %Total number of active excitatory-inhibitory firing rate unit pairs           

dt = 0.0001;        %Time step for simulation
tmax = 10;          %Duration of simulation
t = 0:dt:tmax;      %Simulation time vector
Nt = length(t);

taue = 0.010;       %Time constant for excitatory cells
taui = 0.010;       %Time constant for Inhibitory cells

Wee0 = 1.3;         %Excit.-Excit. connection strength
Wie0 = -0.2;        %Inhib.-Excit. connection strength
Wei0 = 2.9;         %Excit-Excit connection strength
Weix = 0.5;           %Excit.-Inhibitory Cross connection strength
Wiex = -0.1;
Wii0 = -1;          %Inhib-Inhib connection strength

%% Original values
% Wee0 = 1.5;         %Excit.-Excit. connection strength
% Wie0 = -0.4;        %Inhib.-Excit. connection strength
% Wei0 = 2.5;         %Excit-Excit connection strength
% Weix = 0.5;           %Excit.-Inhibitory Cross connection strength
% Wiex = 0;
% Wii0 = -1;  

I0e = 1;            %Excit. applied current
I0i = 17;           %Inhib. applied current

alpha_e = 1;        %Gain of excit. cells
alpha_i = 0.02;     %Gain of inhib. cells
rmax = 100;         %Maximum firing rate

theta_e = 0;%1;    %threshold of activity for e. cells
theta_i = 0;%5;    %threshold of activity for i. cells

%Connection Matricies
Wee = Wee0*eye(N);
Wie = (Wie0-Wiex)*eye(N) + Wiex*ones(N);
Wei = (Wei0-Weix)*eye(N) + Weix*ones(N);
Wii = Wii0*eye(N);

%Current Matricies
Ie = I0e*ones(N,1);
Ii = I0i*ones(N,1);

% Calculate if the parameters fulfills requirments for inhibition stablized
% regime with linear firing rate function.
Wie_tilde = (-Wie0)/(1-Wii0);
I0i_tilde = I0i*Wie_tilde;
off_stable = Wee0 - 1 - Wie_tilde*(Wei0-Weix);
disp(off_stable)

%Rate Matricies
initial_fr = 20;    %Initial firing rate, can be altered.
re = zeros(N,Nt);
ri = zeros(N,Nt);
re(1:M,1) = initial_fr;%50.3;
ri(1:M,1) = 0;

%% noise term, start_t and end_t correspond to the noise term time
% if cut the noise, set start_t and end_t = 0
% if want to add noise to all time, set start_t = 0, end_t = tmax
sigman = 0.005/sqrt(dt);    % initial value is 0.01 
start_t = 6;    % noise start time
end_t = 6.5;    % noise ending time

%% simulation
for i = 2:Nt
    % Remove noise outside of limited range of time (start_t - end_t)
    if (i*dt < start_t || i*dt > end_t)
        Ie_tot = Wee*re(:,i-1) + Wie*ri(:,i-1) + Ie;
        Ii_tot = Wei*re(:,i-1) + Wii*ri(:,i-1) + Ii;
    else
        Ie_tot = Wee*re(:,i-1) + Wie*ri(:,i-1) + Ie + sigman*randn(N,1);
        Ii_tot = Wei*re(:,i-1) + Wii*ri(:,i-1) + Ii + sigman*randn(N,1);        
    end
    Ie_tot = max(Ie_tot,0);
    Ii_tot = max(Ii_tot,0);

    % Calculate excit. firing rate using linear firing rate function;
    % inhib. firing rate using quadratic firing rate function.
    re(:,i) = re(:,i-1) + dt*(-eye(N)*re(:,i-1) + alpha_e*(Ie_tot-theta_e))/taue;   
    ri(:,i) = ri(:,i-1) + dt*(-eye(N)*ri(:,i-1) + alpha_i*(Ii_tot-theta_i).^2).*sign(Ii_tot-theta_i)/taui;
    
    % Firing rate in range of 0-rmax
    re(:,i) = max(re(:,i),0);
    ri(:,i) = max(ri(:,i),0);
    re(:,i) = min(re(:,i),rmax);
    ri(:,i) = min(ri(:,i),rmax);
end

legend_array = zeros(N);

figure(1)
clf
subplot(2,1,1)
plot(t,re), legend("Active", "Non-active"), xlabel("Excit."), ylabel("Firing Rate")
subplot(2,1,2)
plot(t,ri), legend("Active","Non-active"), xlabel("Inhib."), ylabel("Firing Rate")


% add all subfolders of current directory into matlab session search
current_path = pwd;
addpath(genpath(current_path))

disp(is_bistable(N, M, re(:,ceil(0.5/dt:end))))




