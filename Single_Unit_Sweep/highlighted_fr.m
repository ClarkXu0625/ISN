% Coupled Inhibition Stabilized Circuits
%IS regimine with non linear firing rate equations for Excitatory units and
%no cross connections.

% Written for faster simulation
% Single trial that plots the firing rate curve of designated parameters 

clear

%% 
N=1;          %value for how many coupled networks

tmax = 5;
dt = 0.0001;  
tvec = 0:dt:tmax;
rmax = 100;    

frmat_e = zeros(N,numel(tvec));
Imat_e = frmat_e;
frmat_i = zeros(N,numel(tvec));
Imat_i = frmat_i;

%% parameters 
theta_e = 1;    %0:0.1:10;   %threshold of activity for e. cells
theta_i = 5;    %0:0.1:10;   %threshold of activity for i. cells

alpha_e = 0.1; %gain of e. cells
alpha_i = 1;   %gain of i. cells
tao_e = 10e-3; %time constant of e. cells
tao_i = 5e-3;  %time constant of i. cells

%WEI = 3.5;     %connection strength from e. to i. cells 3.5
%WEE = 3;       %connection strength from e. cell to self 3
WII = -6;      %connection strength from i. cell to self
WIE = -3;    %connection strength from i. to e cells
WEIX = 0;     %connection strength from e. cells to i cells from other coupled units. must be changed with N
WEI = 4.9;
WEE = 2.8;
Ii_base=-10;
Ie_base=-1;


%% applied charge
i_stimmat = zeros(N,numel(tvec));%each row is a vector for applied current to each coupled set's inhibitory unit, with # of columns = tmax/dt
Iapp_e = ones(size(i_stimmat))*Ie_base;

%random noise

for i = 1:N
    noisevec = randn(size(tvec))*(dt^(0.5))*15;
    i_stimmat(i,:)= noisevec;
    Iapp_e(i,:) = noisevec + Iapp_e(i,:);
end

%"off switch" 
i_stimmat(1,ceil(2/dt):ceil(2.5/dt))= 15;


Iapp_i = ones(size(i_stimmat))*(Ii_base) + i_stimmat;



% add all subfolders of current directory into matlab session search
current_path = '/Users/apple/Documents/GitHub/Lab/Single_Unit_Sweep';
addpath(genpath(current_path))


%% Simulation, functions are "linI_QE", "linE_QI", and "QEI"
[frmat_i,frmat_e] = ...
    QEI(N, tvec, dt, WEI, WEE, WII, WIE, WEIX, Iapp_i, Iapp_e, theta_i, theta_e, tao_i, tao_e, alpha_e, alpha_i, rmax);

% plot time-firingRate curve for designated paremeter values

clf
figure(2), 
plot(tvec, frmat_e)%, hold on, plot(tvec, frmat_i), legend('e-unit', 'i-unit')
title("WEE==" +num2str(WEE)+" && WEI==" +num2str(WEI)+" && WII==" +num2str(WII))
xlabel("time"), ylabel("firing rate")

disp(std(frmat_e(round(0.2/dt):round(1.9/dt))))
[f, oscillates] = spectrum(frmat_e, tvec, dt);
figure(3), plot(f, oscillates), title("power spectrum of highlighted spot"), xlabel("firing rate (Hz)")
            
        