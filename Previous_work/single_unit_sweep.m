%Coupled Inhibition Stabilized Circuits
%7/26/22
%Code written by Connor Zawacki for the Miller Lab at Brandeis University

%IS regimine with non linear firing rate equations for Excitatory units and
%no cross connections. 

%Looking for an single unit bistability



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
theta_e_vec = 0:0.1:10;   %threshold of activity for e. cells
theta_i_vec = 0:0.1:10;   %threshold of activity for i. cells
alpha_e = 0.1; %gain of e. cells
alpha_i = 1;   %gain of i. cells
tao_e = 10e-3; %time constant of e. cells
tao_i = 5e-3;  %time constant of i. cells

WEI = 3.5;     %connection strength from e. to i. cells 3.5
WEE = 3;       %connection strength from e. cell to self 3
WII = -3;      %connection strength from i. cell to self
WIE = -3.5;    %connection strength from i. to e cells
WEIX = 0;     %connection strength from e. cells to i cells from other coupled units. must be changed with N

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
%}
%"off switch" 
i_stimmat(1,ceil(2/dt):ceil(2.5/dt))= 15;


figure(98), plot(tvec, i_stimmat)


%"on switch" 
frmat_e(1, 1)= 5;
%frmat_i(1:on, 1)= 5;




Iapp_i = ones(size(i_stimmat))*(Ii_base) + i_stimmat;

outputmat = zeros(4,10201);
for i = 1:11
    outputmat(1,(101*(i-1))+1:(101*(i-1))+101)=theta_e_vec(i);
    outputmat(2,(101*(i-1))+1:(101*(i-1))+101)=theta_i_vec;
end
i=0;
for theta_e=theta_e_vec
    for theta_i=theta_i_vec
    i=i+1;

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
        if works==true
            Nss = 1; %follows combination formula, if 1 unit is bistable, any number out of 20 can be active at once
            outputmat(3,i)=outputmat(3,i) + Nss;
            outputmat(4,i)= mean(frmat_e(ceil(0.1/dt):floor(1.99/dt)));
        end

    end
end

figure(99)
imagemat=zeros(101,101);

for i = 1:101
    imagemat(i,:) = outputmat(3,(i-1)*101+1:(i-1)*101+101);
end

x = [0.1,10];
y = [0.1,10];
imagesc(x,y,imagemat);
set(gca,'YDir','normal');
xlabel("Inhibitory Threshold")
ylabel("Excitatory Threshold")
%c = colorbar();
title("Testing theta values for bistability in a single unit")

