%Code written by Connor Zawacki for the Miller Lab at Brandeis University

%Single E-I pair with a nonlinear inhibitory firing rate equation
%exhibiting bistability

%% 
N=1;

tmax = 15;
dt = 0.0001;  
tvec = 0:dt:tmax;
rmax = 100;    %max frate

frmat_e = zeros(N,numel(tvec));
Imat_e = frmat_e;
frmat_i = zeros(N,numel(tvec));
Imat_i = frmat_i;

%% parameters to vary
theta_e = 0;  %threshold of activity for e. cells
theta_i = 0;  %threshold of activity for i. cells
alpha_e = 1;   %gain of e. cells    
alpha_i = 0.02;   %gain of i. cells
tao_e = 0.01;  %time constant of e. cells
tao_i = 0.01;  %time constant of i. cells

WEI = 2.5;     %connection strength from e. to i. cells
WEE = 1.5;       %connection strength from e. cell to self
WII = -1;      %connection strength from i. cell to self
WIE = -0.4;    %connection strength from i. to e cells
WIEX = 0;   %connection strength from e. cells to i cells from other coupled units. must be changed with N


Ii_base=17;
Ie_base=1;

%% applied charge
i_stimmat = zeros(N,numel(tvec));%each row is a vector for applied current to each coupled set's inhibitory unit, with # of columns = tmax/dt


%random noise
for i = 1:N
    noisevec = randn(size(tvec))*(dt^(0.5))*10;
    i_stimmat(i,:)= noisevec;
end


%transition 1 
i_stimmat(1,ceil(5/dt):ceil(5.1/dt))= -20; 


%transition 2 
%insert here

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
figure(3)
for i=1:N
    subplot(1,1,i), plot(tvec,frmat_e(i,:)),  hold on, plot(tvec,frmat_i(i,:)), ylabel("Fr Unit#"+i), xlabel("Time(s)"), hold off;
end
legend("excit.","inhib.")
figure(4)
for i=1:N
    subplot(1,1,i), plot(tvec,Imat_e(i,:)),  hold on, plot(tvec,Imat_i(i,:)),   ylabel("Total input Unit#"+i), xlabel("Time(s)"), hold off;
end
legend("excit.","inhib."),
figure(5)
for i=1:N
    subplot(1,1,i), plot(tvec,Iapp_mat(i,:)), ylabel("Total Iapp Unit#"+i), xlabel("Time(s)"), hold off;
end
