% Coupled Inhibition Stabilized Circuits
%IS regimine with non linear firing rate equations for Excitatory units and
%no cross connections.

% Looking for an single unit bistability

% Exploration on ranges of different parameters.
% This template could be used to simulate various firing rate models (i.e.
% linI-QE, linE-QI, QEI, linEI) by calling matched simulation function file
% Written by Yilin Xu

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
%WII = -3;      %connection strength from i. cell to self
WIE = -3;    %connection strength from i. to e cells
WEIX = 0;     %connection strength from e. cells to i cells from other coupled units. must be changed with N
WEI_vec = 0:0.1:20;
WEE_vec = 0:0.1:20;
WII_vec = -1:-1:-6;
Ii_base=-10;
Ie_base=-1;
Trial = length(WEI_vec);    % number of trials for each parameter

% add all subfolders of current directory into matlab session search
current_path = pwd;
addpath(genpath(current_path))


for m = 1:length(WII_vec)
    WII = WII_vec(m);
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
    
    %"on switch" 
    frmat_e(1, 1)= 5;
    %frmat_i(1:on, 1)= 5;
    
    
    Iapp_i = ones(size(i_stimmat))*(Ii_base) + i_stimmat;
    
    outputmat = zeros(5,Trial,Trial);
    
    for i = 1:Trial
        outputmat(1,:,i) = WEE_vec(i);
        outputmat(2,i,:) = WEI_vec(i);
    end
    
    % Highlighted spot, first place matches the x-axis (WEI), second place 
    % matches the y-axis (WEE).
    highlight = true;   % Whether highlight the spot on figure.
    highlight_spot = [35 41];
    
    
    for i = 1:Trial
        WEE = WEE_vec(i);
        for j = 1:Trial
            WEI = WEI_vec(j);
    
            %% Simulation, functions are "linI_QE", "linE_QI", and "QEI"
            [frmat_i,frmat_e] = ...
                QEI(N, tvec, dt, WEI, WEE, WII, WIE, WEIX, Iapp_i, Iapp_e, theta_i, theta_e, tao_i, tao_e, alpha_e, alpha_i, rmax);
            
            % plot time-firingRate curve for designated paremeter values
            if i==highlight_spot(2) && j==highlight_spot(1)
                figure(96), 
                plot(tvec, frmat_e), title("WEE==" +num2str(WEE_vec(i))+ ...
                    " && WEI=="+num2str(WEI_vec(j))),xlabel("time"), ylabel("firing rate")
            end
    
            Nss = bistability_analysis(frmat_e, tvec, dt);    % calling bistability analysis

            outputmat(3,i,j) = outputmat(3,i,j) + Nss;
            outputmat(4,i,j) = mean(frmat_e(ceil(0.5/dt):floor(1.99/dt)));
            outputmat(5,i,j) = mean(frmat_e(ceil(2.5/dt):floor(4.5/dt)));        
        end
    end
    
    imagemat1 = reshape(outputmat(3,:,:),[Trial,Trial]);
    imagemat2 = reshape(outputmat(4,:,:),[Trial,Trial]);
    imagemat3 = reshape(outputmat(5,:,:),[Trial,Trial]);
    
    % Adding a highlight on the spot that been plotted (used for analysis only)
    if highlight
        imagemat1(highlight_spot(2),highlight_spot(1)) = 3;
    end
    x = [0,20];
    y = [0,20];
    figure(99), subplot(2,3,m), imagesc(x,y,imagemat1);
    set(gca,'YDir','normal'), xlabel("WEI"), ylabel("WEE");
    co = colorbar();
    co.Ticks = [0, 1, 2, 3];
    co.TickLabels = ["single" "bistable" "oscillation" "highlighted"];
    title("Bistability when WII = " + num2str(WII_vec(m)));
    
    
    figure(98), subplot(2,3,m), imagesc(x,y,imagemat2), set(gca,'YDir','normal');
    xlabel("WEI"), ylabel("WEE"), title("E-unit average fr, W_II = " + num2str(WII_vec(m)))
    %c = colorbar;
    %c.Label.String = "firing rate (Hz)";
    
    %figure(97), imagesc(x,y,imagemat3), set(gca,'YDir','normal');
    %xlabel("WEI"), ylabel("WEE"), title("E-unit firing rate after inhibition input given")
    %colorbar
end