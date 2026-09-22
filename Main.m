clc; clear; close all;

% --- Parameters of BPSK Monte Carlo Simulations ---
N = 1e6;                 % Times of Monte Carlo Simulations
SNR_dB = 0:1:15;         % Range of SNR ber bit
N0 = 1;                  % Variance of the Noise

% Parameters of Fading and Random Deployment
K = 5;                   % Rice Factor
h_d = 1;                 % Unit Magnitude Deterministic Complex Scalar
R = 3;                   % Radius of Range of Random Deployment
alpha = 2.2;             % Loss Index

% Initialize 6 arrays to store the simulated bit error rate.
BER_sim_AWGN = zeros(1, length(SNR_dB));
BER_sim_Laplacian = zeros(1, length(SNR_dB));
BER_sim_Rayleigh = zeros(1, length(SNR_dB));
BER_sim_Rician = zeros(1, length(SNR_dB));
BER_sim_Rayleigh_Rand = zeros(1, length(SNR_dB));
BER_sim_Rician_Rand = zeros(1, length(SNR_dB));

% Initialize outage probability related parameters and arrays
target_rate = 1.2;            % C=1.2bps/Hz
gamma_th = 2^target_rate - 1; % outage threshold
gamma_b_linear = 10.^(SNR_dB / 10); 
Outage_sim_Rayleigh = zeros(1, length(SNR_dB));
Outage_sim_Rician = zeros(1, length(SNR_dB));
Outage_sim_Rayleigh_Rand = zeros(1, length(SNR_dB));
Outage_sim_Rician_Rand = zeros(1, length(SNR_dB));


%===========================================================
%                      Simulation Part
%===========================================================

% Signal Sourse
tx_bits = randi([0, 1], 1, N);

fprintf('Simulating...\n');

for k = 1:length(SNR_dB)
    
    gamma_b = 10^(SNR_dB(k) / 10);
    Eb = gamma_b;
    
    % 2. BPSK Modulation
    tx_symbols = (2 * tx_bits - 1) * sqrt(Eb);
    
    % 3. Create the parameters of different channel
    noise_A = AWGN(N, N0);
    noise_L = Laplacian(N, N0); 
    h_s = Rayleigh(N);
    h_R = Rician(N, K, h_d);
    
    % Create Random Distants and Path Loss
    d = R * sqrt(rand(1, N)); % CDF : F(d) = 2d/R^2
    % Path Loss (Amplitude): (d^(-alpha/2))
    path_loss_amp = d.^(-alpha / 2);
    
    % =========================================================
    % Channel 1 & 2: Only with Noise
    % =========================================================
    rx_bits_A = real(tx_symbols + noise_A) > 0;
    BER_sim_AWGN(k) = sum(tx_bits ~= rx_bits_A) / N;
    
    rx_bits_L = (tx_symbols + noise_L) > 0;
    BER_sim_Laplacian(k) = sum(tx_bits ~= rx_bits_L) / N;
    
    % =========================================================
    % Channel 3 & 4: Known distance (d=1) Fading Channels
    % =========================================================
    % Rayleigh
    rx_eq_R = conj(h_s) .* (h_s .* tx_symbols + noise_A);
    rx_bits_R = real(rx_eq_R) > 0;
    BER_sim_Rayleigh(k) = sum(tx_bits ~= rx_bits_R) / N;
    Outage_sim_Rayleigh(k) = sum( gamma_b * abs(h_s).^2 < gamma_th ) / N;

    % Rician
    rx_eq_Ric = conj(h_R) .* (h_R .* tx_symbols + noise_A);
    rx_bits_Ric = real(rx_eq_Ric) > 0;
    BER_sim_Rician(k) = sum(tx_bits ~= rx_bits_Ric) / N;
    Outage_sim_Rician(k) = sum( gamma_b * abs(h_R).^2 < gamma_th ) / N;

    % =========================================================
    % channel 5: Random distance + Rayleigh Fading + AWGN 
    % =========================================================
    % Rayleigh coefficient times path loss amplitude to get combined channel
    h_combined_R = h_s .* path_loss_amp;
    rx_eq_R_rand = conj(h_combined_R) .* (h_combined_R .* tx_symbols + noise_A);
    rx_bits_R_rand = real(rx_eq_R_rand) > 0;
    BER_sim_Rayleigh_Rand(k) = sum(tx_bits ~= rx_bits_R_rand) / N;
    Outage_sim_Rayleigh_Rand(k) = sum( gamma_b * abs(h_combined_R).^2 < gamma_th ) / N;

    % =========================================================
    % channel 6: Random distance + Rician Fading + AWGN 
    % =========================================================
    % Rician coefficient times path loss amplitude to get combined channel
    h_combined_Ric = h_R .* path_loss_amp;
    rx_eq_Ric_rand = conj(h_combined_Ric) .* (h_combined_Ric .* tx_symbols + noise_A);
    rx_bits_Ric_rand = real(rx_eq_Ric_rand) > 0;
    BER_sim_Rician_Rand(k) = sum(tx_bits ~= rx_bits_Ric_rand) / N;
    Outage_sim_Rician_Rand(k) = sum( gamma_b * abs(h_combined_Ric).^2 < gamma_th ) / N;

end
disp('Completed！\n');

% =========================================================
% Theoretical value calculation (all 8)
% =========================================================

% ---------------- Channel 1 & 2 ----------------
BER_theory_AWGN = 0.5 * erfc(sqrt(gamma_b_linear));                                
BER_theory_Laplacian = 0.5 * exp(-sqrt(2 * gamma_b_linear));                       

% Initialize BER and Outage arrays
BER_theory_Rayleigh = zeros(1, length(SNR_dB));
BER_theory_Rician = zeros(1, length(SNR_dB));
BER_theory_Rayleigh_Rand = zeros(1, length(SNR_dB));
BER_theory_Rician_Rand = zeros(1, length(SNR_dB));

Outage_theory_Rayleigh = zeros(1, length(SNR_dB));
Outage_theory_Rician = zeros(1, length(SNR_dB));
Outage_theory_Rayleigh_Rand = zeros(1, length(SNR_dB));
Outage_theory_Rician_Rand = zeros(1, length(SNR_dB));

%---------------- Other Channels -----------------
fprintf('Theoretical value calculation...\n');
for k = 1:length(SNR_dB)
    gb = gamma_b_linear(k);
    
    % --- BER Part ---
    % 1. Rayleigh
    BER_theory_Rayleigh(k) = 0.5 * (1 - sqrt(gb / (1 + gb)));
    
    % 2. Rician
    % use besseli(..., 1) to avoid 0 * Inf = NaN
    pdf_ric = @(x, g) ((K+1)./g) .* exp(-(sqrt((K+1).*x./g) - sqrt(K)).^2) .* besseli(0, 2.*sqrt(K.*(K+1).*x./g), 1);
    % not use Inf, big number like 50*g+50 instead
    fun_ber_ric = @(g) integral(@(x) 0.5.*erfc(sqrt(x)) .* pdf_ric(x, g), 0, 50*g + 50);
    BER_theory_Rician(k) = fun_ber_ric(gb);
    
    % 3. Random Rayleigh
    fun_ber_ray_rand = @(r) (0.5 .* (1 - sqrt((gb.*r.^(-alpha)) ./ (1 + gb.*r.^(-alpha))))) .* (2.*r./R^2);
    BER_theory_Rayleigh_Rand(k) = integral(fun_ber_ray_rand, 1e-6, R);
    
    % 4. Random Rician
    fun_ber_ric_rand = @(r) arrayfun(@(ri) fun_ber_ric(gb.*ri.^(-alpha)) .* (2.*ri./R^2), r);
    BER_theory_Rician_Rand(k) = integral(fun_ber_ric_rand, 1e-6, R);
    
    % --- Outage ---
    % 1. Fixed Rayleigh
    Outage_theory_Rayleigh(k) = 1 - exp(-gamma_th / gb);
    
    % 2. Fixed Rician
    Outage_theory_Rician(k) = 1 - marcumq(sqrt(2*K), sqrt(2*(K+1)*gamma_th / gb));
    
    % 3. Random Rayleigh 
    fun_outage_ray_rand = @(r) (1 - exp(-gamma_th ./ (gb.*r.^(-alpha)))) .* (2.*r./R^2);
    Outage_theory_Rayleigh_Rand(k) = integral(fun_outage_ray_rand, 1e-6, R);
    
    % 4. Random Rician 
    fun_outage_ric_rand = @(r) (1 - marcumq(sqrt(2*K), sqrt(2*(K+1)*gamma_th ./ (gb.*r.^(-alpha))))) .* (2.*r./R^2);
    Outage_theory_Rician_Rand(k) = integral(fun_outage_ric_rand, 1e-6, R);
end


% =========================================================
% Plot Part
% =========================================================

% === Figure 1: Only Noise BER ===
figure(1);
% Plot theortical value (curve)
semilogy(SNR_dB, BER_theory_AWGN, 'b-', 'LineWidth', 1.5); hold on;
semilogy(SNR_dB, BER_theory_Laplacian, 'r-', 'LineWidth', 1.5);
% PLot simulation value (point)
semilogy(SNR_dB, BER_sim_AWGN, 'bo', 'MarkerSize', 6, 'LineWidth', 1.5); 
semilogy(SNR_dB, BER_sim_Laplacian, 'r*', 'MarkerSize', 6, 'LineWidth', 1.5); 
grid on;
xlabel('SNR (\gamma_b) [dB]', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('Plot 1: BER in Noise-Only Environments', 'FontSize', 14);
legend('Analytical AWGN', 'Analytical Laplacian', 'Simulated AWGN', 'Simulated Laplacian', 'Location', 'southwest');
axis([0 15 1e-6 1]); 

% === Figure 2: Fading channel BER with known & random distance ===
figure(2);
% Plot theortical value (curve)
semilogy(SNR_dB, BER_theory_Rayleigh, 'k-', 'LineWidth', 1.5); hold on;
semilogy(SNR_dB, BER_theory_Rician, 'm-', 'LineWidth', 1.5); 
semilogy(SNR_dB, BER_theory_Rayleigh_Rand, 'g-', 'LineWidth', 1.5); 
semilogy(SNR_dB, BER_theory_Rician_Rand, 'c-', 'LineWidth', 1.5); 
% PLot simulation value (point)
semilogy(SNR_dB, BER_sim_Rayleigh, 'ks', 'MarkerSize', 6, 'LineWidth', 1.5); 
semilogy(SNR_dB, BER_sim_Rician, 'md', 'MarkerSize', 6, 'LineWidth', 1.5); 
semilogy(SNR_dB, BER_sim_Rayleigh_Rand, 'g^', 'MarkerSize', 6, 'LineWidth', 1.5); 
semilogy(SNR_dB, BER_sim_Rician_Rand, 'cv', 'MarkerSize', 6, 'LineWidth', 1.5); 
grid on;
xlabel('SNR (\gamma_b) [dB]', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('Plot 2: BER in Fading & Random Deployment', 'FontSize', 14);
legend('Theory Fixed Rayleigh', 'Theory Fixed Rician', 'Theory Rand Rayleigh', 'Theory Rand Rician', ...
       'Sim Fixed Rayleigh', 'Sim Fixed Rician', 'Sim Rand Rayleigh', 'Sim Rand Rician', 'Location', 'southwest');
axis([0 15 1e-4 1]);

% === Figure 3: Outage probability in 4 types of channel ===
figure(3);
% Plot theortical value (curve)
semilogy(SNR_dB, Outage_theory_Rayleigh, 'k-', 'LineWidth', 1.5); hold on;
semilogy(SNR_dB, Outage_theory_Rician, 'm-', 'LineWidth', 1.5); 
semilogy(SNR_dB, Outage_theory_Rayleigh_Rand, 'g-', 'LineWidth', 1.5); 
semilogy(SNR_dB, Outage_theory_Rician_Rand, 'c-', 'LineWidth', 1.5); 
% PLot simulation value (point)
semilogy(SNR_dB, Outage_sim_Rayleigh, 'ks', 'MarkerSize', 6, 'LineWidth', 1.5); 
semilogy(SNR_dB, Outage_sim_Rician, 'md', 'MarkerSize', 6, 'LineWidth', 1.5); 
semilogy(SNR_dB, Outage_sim_Rayleigh_Rand, 'g^', 'MarkerSize', 6, 'LineWidth', 1.5); 
semilogy(SNR_dB, Outage_sim_Rician_Rand, 'cv', 'MarkerSize', 6, 'LineWidth', 1.5); 
grid on;
xlabel('Average SNR (\gamma_b) [dB]', 'FontSize', 12);
ylabel('Outage Probability (P_{out})', 'FontSize', 12);
title('Plot 3: Outage Probability in Fading Environments', 'FontSize', 14);
legend('Theory Fixed Rayleigh', 'Theory Fixed Rician', 'Theory Rand Rayleigh', 'Theory Rand Rician', ...
       'Sim Fixed Rayleigh', 'Sim Fixed Rician', 'Sim Rand Rayleigh', 'Sim Rand Rician', 'Location', 'southwest');
axis([0 15 1e-3 1]);


% =========================================================
% Laplacian Noise PDF & CDF Therotical vs Simulation
% =========================================================

figure ('Name', 'Laplacian Verification', 'Position', [150, 150, 1000, 400]);

scale_b = sqrt(N0 / 2);         % Scale Factor
x_val = linspace(-5, 5, 1000);  % x axis

% --- Therotical calculation ---
pdf_theory = (1 / (2 * scale_b)) * exp(-abs(x_val) / scale_b);
cdf_theory = zeros(size(x_val));
cdf_theory(x_val <= 0) = 0.5 * exp(x_val(x_val <= 0) / scale_b);
cdf_theory(x_val > 0)  = 1 - 0.5 * exp(-x_val(x_val > 0) / scale_b);

noise_samples = real(noise_L); % only real part

%-------- PDF ----------
subplot(1, 2, 1);
% simulation
histogram(noise_samples, 1000, 'Normalization', 'pdf', 'FaceColor', [0.7 0.7 0.9], 'EdgeColor', 'none');
hold on;
% Therotical
plot(x_val, pdf_theory, 'r-', 'LineWidth', 2);
title('Laplacian PDF: Simulation vs Theory');
xlabel('Noise Amplitude (x)');
ylabel('Probability Density f(x)');
legend('Simulated Histogram', 'Theoretical Curve');
xlim([-5 5]); grid on;

% ------  CDF -------
subplot(1, 2, 2);
% simulation
[f_emp, x_emp] = ecdf(noise_samples);
plot(x_emp(1:1000:end), f_emp(1:1000:end), 'ko', 'MarkerSize', 5); 
hold on;
% therotical
plot(x_val, cdf_theory, 'r-', 'LineWidth', 2);
title('Laplacian CDF: Simulation vs Theory');
xlabel('Noise Amplitude (x)');
ylabel('Cumulative Probability F(x)');
legend('Simulated (Empirical CDF)', 'Theoretical Curve', 'Location', 'southeast');
xlim([-5 5]); grid on;


% ---------------------------------------------------------
% Command Window Output
% ---------------------------------------------------------

% BER Form
fprintf('\n==========================================================================================================\n');
fprintf('                                        BER Simulation Value \n');
fprintf('==========================================================================================================\n');
fprintf(' SNR(dB) |     AWGN     |  Laplacian   | Fixed Rayleigh | Fixed Rician | Rand Rayleigh | Rand Rician \n');
fprintf('----------------------------------------------------------------------------------------------------------\n');
for k = 1:length(SNR_dB)
    fprintf('   %2d    |  %1.4e  |  %1.4e  |   %1.4e   |  %1.4e  |   %1.4e  |  %1.4e  \n', ...
        SNR_dB(k), BER_sim_AWGN(k), BER_sim_Laplacian(k), BER_sim_Rayleigh(k), BER_sim_Rician(k), BER_sim_Rayleigh_Rand(k), BER_sim_Rician_Rand(k));
end
fprintf('==========================================================================================================\n');

% Outage Probability Form
fprintf('\n==================================================================================\n');
fprintf('                           Outage Probability Simulation Value \n');
fprintf('==================================================================================\n');
fprintf(' SNR(dB) | Fixed Rayleigh | Fixed Rician | Rand Rayleigh | Rand Rician \n');
fprintf('----------------------------------------------------------------------------------\n');
for k = 1:length(SNR_dB)
    fprintf('   %2d    |   %1.4e   |  %1.4e  |   %1.4e  |  %1.4e  \n', ...
        SNR_dB(k), Outage_sim_Rayleigh(k), Outage_sim_Rician(k), Outage_sim_Rayleigh_Rand(k), Outage_sim_Rician_Rand(k));
end
fprintf('==================================================================================\n\n');


% =========================================================================
%                              Percentage Error
% =========================================================================
disp(' ');
disp('=========================================================================================================');
disp('                              Percentage Error Analysis (Simulation vs Theory)                           ');
disp('=========================================================================================================');

% Calculate err in all 6 channels and all SNR
% Plus eps to prevent NaN
err_BER_AWGN       = abs(BER_sim_AWGN - BER_theory_AWGN) ./ (BER_theory_AWGN + eps) * 100;
err_BER_Laplacian  = abs(BER_sim_Laplacian - BER_theory_Laplacian) ./ (BER_theory_Laplacian + eps) * 100;
err_BER_Rayleigh   = abs(BER_sim_Rayleigh - BER_theory_Rayleigh) ./ (BER_theory_Rayleigh + eps) * 100;
err_BER_Rician     = abs(BER_sim_Rician - BER_theory_Rician) ./ (BER_theory_Rician + eps) * 100;
err_BER_Ray_Rand   = abs(BER_sim_Rayleigh_Rand - BER_theory_Rayleigh_Rand) ./ (BER_theory_Rayleigh_Rand + eps) * 100;
err_BER_Ric_Rand   = abs(BER_sim_Rician_Rand - BER_theory_Rician_Rand) ./ (BER_theory_Rician_Rand + eps) * 100;

err_Out_Rayleigh   = abs(Outage_sim_Rayleigh - Outage_theory_Rayleigh) ./ (Outage_theory_Rayleigh + eps) * 100;
err_Out_Rician     = abs(Outage_sim_Rician - Outage_theory_Rician) ./ (Outage_theory_Rician + eps) * 100;
err_Out_Ray_Rand   = abs(Outage_sim_Rayleigh_Rand - Outage_theory_Rayleigh_Rand) ./ (Outage_theory_Rayleigh_Rand + eps) * 100;
err_Out_Ric_Rand   = abs(Outage_sim_Rician_Rand - Outage_theory_Rician_Rand) ./ (Outage_theory_Rician_Rand + eps) * 100;

% printf form
fprintf(' SNR |     BER Percentage Error (%%)     |     Outage Percentage Error (%%)    \n');
fprintf('(dB) | AWGN  Lapl  Rayl  Rici  Ray(R) Ric(R) |  Rayl  Rici  Ray(R) Ric(R) \n');
disp('---------------------------------------------------------------------------------------------------------');

% all value
for k = 1:length(SNR_dB)
    fprintf(' %2d  | %5.2f %5.2f %5.2f %5.2f  %5.2f  %5.2f | %5.2f %5.2f  %5.2f  %5.2f \n', ...
        SNR_dB(k), ...
        err_BER_AWGN(k), err_BER_Laplacian(k), err_BER_Rayleigh(k), err_BER_Rician(k), err_BER_Ray_Rand(k), err_BER_Ric_Rand(k), ...
        err_Out_Rayleigh(k), err_Out_Rician(k), err_Out_Ray_Rand(k), err_Out_Ric_Rand(k));
end
disp('=========================================================================================================');

% mean of values
fprintf(' MEAN| %5.2f %5.2f %5.2f %5.2f  %5.2f  %5.2f | %5.2f %5.2f  %5.2f  %5.2f \n', ...
    mean(err_BER_AWGN, 'omitnan'), mean(err_BER_Laplacian, 'omitnan'), mean(err_BER_Rayleigh, 'omitnan'), ...
    mean(err_BER_Rician, 'omitnan'), mean(err_BER_Ray_Rand, 'omitnan'), mean(err_BER_Ric_Rand, 'omitnan'), ...
    mean(err_Out_Rayleigh, 'omitnan'), mean(err_Out_Rician, 'omitnan'), mean(err_Out_Ray_Rand, 'omitnan'), ...
    mean(err_Out_Ric_Rand, 'omitnan'));
disp('=========================================================================================================');
