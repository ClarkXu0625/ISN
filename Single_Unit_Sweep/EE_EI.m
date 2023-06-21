%Coupled Inhibition Stabilized Circuits

%IS regimine with non linear firing rate equations for Excitatory units and
%no cross connections. 

% Looking for an single unit bistability

% exploration on ranges of different parameters.


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
WII = -3;      %connection strength from i. cell to self
WIE = -3.5;    %connection strength from i. to e cells
WEIX = 0;     %connection strength from e. cells to i cells from other coupled units. must be changed with N
WEI_vec = 0:0.1:10;
WEE_vec = 0:0.1:10;
Ii_base=-10;
Ie_base=-1;
Trial = length(WEI_vec);    % number of trials for each parameter

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

outputmat = zeros(4,Trial,Trial);
for i = 1:11    
    outputmat(1,(101*(i-1))+1:(101*(i-1))+101)=WEE_vec(i);
    outputmat(2,(101*(i-1))+1:(101*(i-1))+101)=WEI_vec;
end

for i = 1:Trial
    outputmat(1,:,i) = WEE_vec(i);
    outputmat(2,i,:) = WEI_vec(i);
end

% Highlighted spot, first place matches the x-axis (WEI), second place 
% matches the y-axis (WEE).
highlight = [36 33];
for i = 1:Trial
    WEE = WEE_vec(i);
    for j = 1:Trial
        WEI = WEI_vec(j);
    
        %% simulation
        for t = 2:numel(tvec)
                Imat_i(:,t) = WEI*frmat_e(:,t-1) + WII*frmat_i(:,t-1) + Iapp_i(:,t) + WEIX*(sum(frmat_e(:,t-1)) - frmat_e(:,t-1));
                Imat_e(:,t) = WEE*frmat_e(:,t-1) + WIE*frmat_i(:, t-1) + Iapp_e(:,t);
        
                frmat_e(:,t) = frmat_e(:,t-1) + (dt/tao_e).*(-frmat_e(:,t-1) + alpha_e*(Imat_e(:,t)-theta_e).^2 .*sign(Imat_e(:,t)-theta_e));
                frmat_i(:,t) = frmat_i(:,t-1) + (dt/tao_i).*(-frmat_i(:,t-1) + alpha_i*(Imat_i(:,t)-theta_i));
                
                frmat_i(:,t)=min(frmat_i(:,t),rmax);
                frmat_e(:,t)=min(frmat_e(:,t),rmax);
                frmat_i(:,t)=max(frmat_i(:,t),0);
                frmat_e(:,t)=max(frmat_e(:,t),0);
        
        end
        
        % plot time-firingRate curve for designated paremeter values
        if i==highlight(2) && j==highlight(1)
            figure(97), 
            plot(tvec, frmat_e), title("WEE==" +num2str(WEE_vec(i))+ ...
                " && WEI=="+num2str(WEI_vec(j))),xlabel("time"), ylabel("firing rate")
        end
        works = true;       
        if frmat_e(ceil(1.5/dt))<0.2
            works=false;
        end
        if frmat_e(ceil(1.5/dt))>90
            works=false;
        end
        if frmat_e(ceil(4/dt))>0.2
            works=false;
        end
        % Test whether there in intrinsic oscillation
        if std(frmat_e(ceil(0.1/dt):ceil(1.9/dt)))>1
            works=false;
        end

        if works==true
            Nss = 1; %follows combination formula, if 1 unit is bistable, any number out of 20 can be active at once
            outputmat(3,i,j)=outputmat(3,i,j) + Nss;
            outputmat(4,i,j)= mean(frmat_e(ceil(0.1/dt):floor(1.99/dt)));
        end
    
    end
end

imagemat1 = reshape(outputmat(3,:,:),[Trial,Trial]);
imagemat2 = reshape(outputmat(4,:,:),[Trial,Trial]);

% Adding a highlight on the spot that been plotted to see if the right
% value been analyzed
imagemat1(highlight(2),highlight(1)) = imagemat1(highlight(2),highlight(1)) + 2;

x = [0.1,10];
y = [0.1,10];
figure(99), imagesc(x,y,imagemat1);
set(gca,'YDir','normal'), xlabel("WEI"), ylabel("WEE");
%c = colorbar();
title("Testing theta values for bistability in a single unit")


figure(98), imagesc(x,y,imagemat2), set(gca,'YDir','normal');
xlabel("WEI"), ylabel("WEE");
colorbar