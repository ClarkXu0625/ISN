function [re, ri] = linE_QI(N, M, Nt, dt, Wee0, Wie0, Wiex, Wii0, Wei0, Weix, sigman, Ie, Ii, taue, taui, alpha_e, alpha_i)

    %Connection Matricies
    Wee = Wee0*eye(N);
    Wie = (Wie0-Wiex)*eye(N) + Wiex*ones(N);
    Wei = (Wei0-Weix)*eye(N) + Weix*ones(N);
    Wii = Wii0*eye(N);

    %Rate Matricies
    re = zeros(N,Nt);
    ri = zeros(N,Nt);
    re(1:M,1) = 15;
    ri(1:M,1) = 0;

    for i = 2:Nt
        Ie_tot = (Wee*re(:,i-1)+ Wie*ri(:,i-1) + Ie + sigman*randn(N,1));
        Ie_tot = max(Ie_tot,0);
        re(:,i) = re(:,i-1) + dt*(-eye(N)*re(:,i-1) + alpha_e*Ie_tot)/taue;
        Ii_tot = (Wei*re(:,i-1) + Wii*ri(:,i-1) + Ii + sigman*randn(N,1));
        Ii_tot = max(Ii_tot,0);
        ri(:,i) = ri(:,i-1) + dt*(-eye(N)*ri(:,i-1)+alpha_i*Ii_tot.^2 )/taui;
        
        re(:,i) = max(re(:,i),0);
        ri(:,i) = max(ri(:,i),0);
    end