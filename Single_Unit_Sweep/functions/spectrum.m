function [f, power_spectrum, oscillates] = spectrum(fr, t, dt)
    % Cut the firing rate within 0-2 second
    fr = fr(1:round(2/dt));
    t = t(1:round(2/dt));

    f = 0:1:100;    % firing rate  
    A = zeros(size(f));
    B = zeros(size(f));

    %% Fourier Transform
    for i = 1:length(f)
        A(i) = mean(sin(2*pi*f(i)*t).*fr);
        B(i) = mean(cos(2*pi*f(i)*t).*fr);
    end
    power_spectrum = A.^2 + B.^2;
    power_spectrum = power_spectrum(2:end);
    f = f(2:end);

    %% Detect whether the given firing rate contains oscillations
    oscillates = false;
    if (max(power_spectrum)>0.05)
        oscillates = true;
    end
end