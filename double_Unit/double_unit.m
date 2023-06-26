% Network with N total inhibition-stabilized excitatory-inhibitory units. 
% M of those pairs we attempt to make "active", and N-M pairs we try to 
% suppress.

%In particular, this simulation displays bistability within a single 
%unit given a quadratic firing rate equation for the inhibitory units

%No WEIX cross connection
%No WIEX cross connection
%nonlinear inhibitory unit Firing rate
%linear excitatory unit Firing rate

clear
N = 2;              %Total excitatory-inhibitory firing rate unit pairs
M = 1;              %Total number of active excitatory-inhibitory firing rate unit pairs           

dt = 0.0001;        %Time step for simulation
tmax = 10;          %Duration of simulation
t = 0:dt:tmax;      %Simulation time vector
Nt = length(t);

taue = 0.010;       %Time constant for excitatory cells
taui = 0.010;       %Time constant for Inhibitory cells

Wee0 = 1.5;         %Excit.-Excit. connection strength
Wie0 = -0.4;        %Inhib.-Excit. connection strength
Wei0 = 2.5;         %Excit-Excit connection strength
Weix = 0.5;           %Excit.-Inhibitory Cross connection strength
Wiex = 0;
Wii0 = -1;          %Inhib-Inhib connection strength

I0e = 1;            %Excit. applied current
I0i = 17;           %Inhib. applied current

alpha_e = 1;        %Gain of excit. cells
alpha_i = 0.02;     %Gain of inhib. cells

%Connection Matricies
Wee = Wee0*eye(N);
Wie = (Wie0-Wiex)*eye(N) + Wiex*ones(N);
Wei = (Wei0-Weix)*eye(N) + Weix*ones(N);
Wii = Wii0*eye(N);

%Current Matricies
Ie = I0e*ones(N,1);
Ii = I0i*ones(N,1);

Wie_tilde = (-Wie0)/(1-Wii0);
I0i_tilde = I0i*Wie_tilde;
off_stable = Wee0 - 1 - Wie_tilde*(Wei0-Weix);
disp(off_stable)

%Rate Matricies
re = zeros(N,Nt);
ri = zeros(N,Nt);
re(1:M,1) = 15;
ri(1:M,1) = 0;

sigman = 0.01/sqrt(dt);

%simulation
for i = 2:Nt
    Ie_tot = (Wee*re(:,i-1)+ Wie*ri(:,i-1) + Ie + sigman*randn(N,1));
    Ie_tot = max(Ie_tot,0);
    re(:,i) = re(:,i-1) + dt*(-eye(N)*re(:,i-1) + alpha_e*Ie_tot )/taue;
    Ii_tot = (Wei*re(:,i-1) + Wii*ri(:,i-1) + Ii + sigman*randn(N,1));
    Ii_tot = max(Ii_tot,0);
    ri(:,i) = ri(:,i-1) + dt*(-eye(N)*ri(:,i-1)+alpha_i*Ii_tot.^2 )/taui;
    
    re(:,i) = max(re(:,i),0);
    ri(:,i) = max(ri(:,i),0);
end

figure(1)
clf
subplot(2,1,1)
plot(t,re), legend("Active","Non-active"), xlabel("Excit."), ylabel("Firing Rate")
subplot(2,1,2)
plot(t,ri), legend("Active","Non-active"), xlabel("Inhib."), ylabel("Firing Rate")

disp(is_bistable(N, re(:,ceil(1/dt:end)),ri(:,ceil(1/dt:end))))


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

