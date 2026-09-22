function n_Laplace = Laplacian(N, N0)
%   N  - number of sample needed
%   N0 - average noise power(variance)

    b = sqrt(N0 / 2);
    U = rand(1, N);
    V = U - 0.5;
    
    % use inverse CDF funtion to describe laplacian noise
    n_Laplace = -b * sign(V) .* log(1 - 2 * abs(V));
    
end