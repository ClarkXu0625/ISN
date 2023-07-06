function Nss = bistability_analysis(fr, t, dt)
% Function bistability_analysis analyzes the firing rate curve. Nss value
% is 1 when there's bistability, and 0 when the unit is quiescense or
% always firing at high firing rate. Nss returns 2 when they show
% oscillatory activity. It could be confusing to distinguish noise and
% oscillations simply analyzing the standard deviation of firing rate, so
% it would also call spectrum function when the standard deviation is above
% certain levels (in this case 0.1) to save running time
    
    Nss = 1;
    if fr(ceil(1.5/dt))<0.2
        Nss = 0;
    end
    if fr(ceil(1.5/dt))>90
        Nss = 0;
    end
    if fr(ceil(4/dt))>0.2
        Nss = 0;
    end
    % Test whether there in intrinsic oscillation
    % Stops at 0.5 second when neurons are mostly settled down
    if std(fr(ceil(0.5/dt):ceil(1.9/dt)))>0.1
        [~, ~, oscillates] = spectrum(fr, t, dt);
        if oscillates
            Nss = 2;
        end
    end

end