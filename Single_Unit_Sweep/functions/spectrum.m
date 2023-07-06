function [f, power_spectrum, oscillates] = spectrum(fr, t)
    f = 0:1:100;
    A = zeros(size(f));
    B = zeros(size(f));
    for i = 1:length(f)
        A(i) = mean(sin(2*pi*f(i)*t).*fr);
        B(i) = mean(cos(2*pi*f(i)*t).*fr);
    end
    power_spectrum = A.^2 + B.^2;
    power_spectrum = power_spectrum(2:end);
    f = f(2:end);

    oscillates = true;

end