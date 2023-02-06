% This code is to written to familiarize with firing-rate model and testing
% for set of parameters
% Data is from Spatiotemporal Discrimination in Attractor Networks with
% Short-term Synaptic Plasticity
%%%%%%%
N_E=100;    % population of excitatory cells
N_I=1;      % population of inhibitory cells

% Time constants
tau_r = 0.100;
tau_sE = 0.050;
tau_sI = 0.005;
tau_D = 0.500;

delta_E = 1;    % sensitivity current
delta_I = 3;
thresholdE = 6;

Inpu