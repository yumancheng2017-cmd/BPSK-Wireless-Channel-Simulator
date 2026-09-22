function n_AWGN = AWGN(N, N0)
%   N  - number of sample needed
%   N0 - average noise power(variance)

    % separate variance for real and imaginary part
    scale = sqrt(N0 / 2);
    
    n_re = scale * randn(1, N);
    n_im = scale * randn(1, N);
    
    n_AWGN = n_re + 1i * n_im;
    
end