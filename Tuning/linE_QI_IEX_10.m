%Coupled Inhibition Stabilized Circuits
%6/28/22
%Code written by Connor Zawacki for the Miller Lab at Brandeis University

%IS regimine with non linear firing rate equations for Inhibitory units and
%I-E cross connections. 

%min number of active pairs: 3
%max number of active pairs: 5
%investigate stability! 
% 


%% 
N=10;          %value for how many coupled networks
on = 6;
tmax = 15;
dt = 0.1e-3;  
tvec = 0:dt:tmax;
rmax = 100;    %max frate

frmat_e = zeros(N,numel(tvec));
Imat_e = frmat_e;
frmat_i = zeros(N,numel(tvec));
Imat_i = frmat_i;

%% parameters to vary
theta_e = -5;  %threshold of activity for e. cells
theta_i = 10;   %threshold of activity for i. cells
alpha_e = 1;   %gain of e. cells
alpha_i = 1;   %gain of i. cells
tao_e = 5e-3;  %time constant of e. cells
tao_i = 5e-3;  %time constant of i. cells

WEI = 2.5;     %connection strength from e. to i. cells
WEE = 2;       %connection strength from e. cell to self
WII = -2;      %connection strength from i. cell to self
WIE = -1.5;    %connection strength from i. to e cells
WIEX = -0.607;     %connection strength from e. cells to i cells from other coupled units. must be changed with N

Ii_base=20;
Ie_base=25;

%% applied charge
i_stimmat = zeros(N,numel(tvec));%each row is a vector for applied current to each coupled set's inhibitory unit, with # of columns = tmax/dt

%random noise
for i = 1:N
    noisevec = randn(size(tvec))*(dt^(0.5))*5;
    i_stimmat(i,:)= noisevec;
end

%transition 1
i_stimmat(1:on,ceil(1/dt):ceil(1.1/dt))= -4;    


Iapp_e = ones(size(tvec))*Ie_base;
Iapp_mat = ones(size(i_stimmat))*(Ii_base) + i_stimmat;

%% simulation
for t = 2:numel(tvec)
    for c = 1:N
        Imat_i(c,t) = WEI*frmat_e(c,t-1) + WII*frmat_i(c,t-1) + Iapp_mat(c,t);
        Imat_e(c,t) = WEE*frmat_e(c,t-1) + WIE*frmat_i(c, t-1) + WIEX*(sum(frmat_i(:,t-1)) - frmat_i(c,t-1)) + Iapp_e(t);

        frmat_e(c,t) = frmat_e(c,t-1) + (dt/tao_e)*(-frmat_e(c,t-1) + alpha_e*(Imat_e(c,t)-theta_e));
        frmat_i(c,t) = frmat_i(c,t-1) + (dt/tao_i)*(-frmat_i(c,t-1) + alpha_i*(Imat_i(c,t)-theta_i)^2 *sign(Imat_i(c,t)-theta_i));

        if frmat_i(c,t)>rmax
            frmat_i(c,t)=rmax;
        end
        if frmat_i(c,t)<0
            frmat_i(c,t)=0;
        end
        if frmat_e(c,t)>rmax
            frmat_e(c,t)=rmax;
        end
        if frmat_e(c,t)<0
            frmat_e(c,t)=0;
        end
    end
end

%% preliminary figures
figure(1)
subplot(2,1,1), plot(tvec, frmat_e), title("E unit");
subplot(2,1,2), plot(tvec, frmat_i), title("I unit");

%figure(2)
%subplot(2,1,1), plot(tvec, Imat_e), title("E unit ");
%subplot(2,1,2), plot(tvec, Imat_i), title("I unit ");

%figure(3)
for i=1:N
    %subplot(4,3,i), plot(tvec,frmat_e(i,:)),  hold on, plot(tvec,frmat_i(i,:)), ylabel("Fr Unit#"+i), xlabel("Time(s)"), hold off;
end
%legend("excit.","inhib.")
%figure(4)
for i=1:N
    %subplot(4,3,i), plot(tvec,Imat_e(i,:)),  hold on, plot(tvec,Imat_i(i,:)),   ylabel("Total input Unit#"+i), xlabel("Time(s)"), hold off;
end
%legend("excit.","inhib."),
%figure(5)
for i=1:N
    %subplot(4,3,i), plot(tvec,Iapp_mat(i,:)), ylabel("Total Iapp Unit#"+i), xlabel("Time(s)"), hold off;
end
