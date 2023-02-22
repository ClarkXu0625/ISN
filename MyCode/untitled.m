% Network with N total inhibition-stabilized excitatory-inhibitory units. 
% M of those pairs we attempt to make "active", and N-M pairs we try to 
% suppress.

%Originally Written by Paul Miller
%Modified and Commented by Connor Zawacki

%WEIX cross connection
%No WIEX cross connection
%linear inhibitory unit Firing rate
%nonlinear excitatory unit Firing rate

clear
N = 10;         %Total excitatory-inhibitory firing rate unit pairs
M = 5;         %Total number of active excitatory-inhibitory firing rate unit pairs

dt = 0.0001;    %Time step for simulation
tmax = 6;       %Duration of Simulation
t = 0:dt:tmax;  %Simulation time vector
Nt = length(t);

taue = 0.010;   %Time constant for excitatory cells
taui = 0.005;   %Time constant for inhibitory cells

Wee0 = 3;       %Excit.-Excit. connection strength
Wie0 = -3;      %Inhib.-Excit. connection strength
Wei0 = 3.5;     %Excit-Inhib connection strength
Weix = 0.025;   %Excit.-Inhibitory Cross connection strength
Wiex = 0.0;     %Inhibitory.-Excitatory Cross connection strength
Wii0 = -3;      %Inhib-Inhib connection strength

I0e = -1;       %Excit. applied current
I0i = -10;      %Inhib. applied current

reh = 30;       %Initial frate of excit cells50
rih = 45;       %Initial frate of inhib cells

alpha_e = 0.1;  %Gain of excit. cells
alpha_i = 1;    %Gain of inhib. cells

%Connection matricies
Wee = Wee0*eye(N);
Wie = Wie0*eye(N);
Wie = (Wie0-Wiex)*eye(N) + Wiex*ones(N);
Wei = (Wei0-Weix)*eye(N) + Weix*ones(N);
Wii = Wii0*eye(N);

%Current Matricies
Ie = I0e*ones(N,1);
Ii = I0i*ones(N,1);

Wei_tilde = Wei0/(1-Wii0);
Wie_tilde = (-Wie0)/(1-Wii0);
Wiex_tilde = (-Wiex)/(1-Wii0);
I0i_tilde = I0i*Wie_tilde;

%Rate Matricies
re = zeros(N,Nt);   
ri = zeros(N,Nt);   
re(1:M,1) = reh;
ri(1:M,1) = rih;

sigman = 0.002/sqrt(dt);

%simulation
for i = 2:Nt
    Ie_tot = (Wee*re(:,i-1)+ Wie*ri(:,i-1) + Ie + sigman*randn(N,1));
    Ie_tot = max(Ie_tot,0);
    re(:,i) = re(:,i-1) + dt*(-eye(N)*re(:,i-1) + alpha_e*Ie_tot.^2 )/taue;
    Ii_tot = (Wei*re(:,i-1) + Wii*ri(:,i-1) + Ii + sigman*randn(N,1));
    Ii_tot = max(Ii_tot,0);
    ri(:,i) = ri(:,i-1) + dt*(-eye(N)*ri(:,i-1)+alpha_i*Ii_tot )/taui;
    re(:,i) = max(re(:,i),0);
    ri(:,i) = max(ri(:,i),0);
end

figure(1)
clf
subplot(2,1,1)
plot(t,re), ylabel("E. Firing_rates")
subplot(2,1,2)
plot(t,ri), ylabel("I. Firing_rates"), xlabel("Time")

