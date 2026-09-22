function h_s = Rayleigh(N)
%   N  - number of sample needed
%   coefficient of channel in Rayleigh fading
    % average value = 0，variance = 1 (real and imaginary each has 0.5)
    h_s = sqrt(0.5) * (randn(1, N) + 1i * randn(1, N));
end