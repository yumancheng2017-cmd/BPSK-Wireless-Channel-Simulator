function h_R = Rician(N, K, h_d)
%   N   - number of sample needed
%   K   - Rician factor (linear)
%   h_d - Direct path component
   
    %   Create scattered components (zero-mean Gaussian)
    h_s = sqrt(0.5) * (randn(1, N) + 1i * randn(1, N));
    
    h_R = sqrt(K / (K + 1)) * h_d + sqrt(1 / (K + 1)) * h_s;
    
end