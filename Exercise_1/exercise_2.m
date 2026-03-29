%% =========================================================
%  Simulation and Modelling of Dynamic Systems
%  Lab 01 - Exercise 2
%% =========================================================

clear; clc; close all;

%% Load data from Exercise 1
% Loads: t, x1, x2, u_out, m, k, c, dt
load('simulation_data.mat');
 
fprintf('===== Data loaded from simulation_data.mat =====\n');
fprintf('True parameters: m=%.4f, k=%.4f, c=%.4f\n\n', m, k, c);

%% =========================================================
%  Exercise 2a — Least Squares Estimation (Ts = 0.05 s)
%% =========================================================
fprintf('===== Exercise 2a =====\n');
Ts = 0.05;  % sampling period [s]
 
% Sample the simulation data at every Ts seconds
step = round(Ts / dt);
idx  = 1 : step : length(t);
x1_samp = x1(idx);
x2_samp = x2(idx);
u_samp  = u_out(idx);
 
% Estimate acceleration via numerical differentiation
x2_dot = gradient(x2_samp, Ts);
 
% Apply Least Squares
[m_hat, c_hat, k_hat] = least_squares(x1_samp, x2_samp, x2_dot, u_samp);
 
fprintf('True parameters:  m=%.4f, c=%.4f, k=%.4f\n', m, c, k);
fprintf('Estimated params: m=%.4f, c=%.4f, k=%.4f\n', m_hat, c_hat, k_hat);
fprintf('Errors:           m=%.4f, c=%.4f, k=%.4f\n\n', ...
        abs(m-m_hat), abs(c-c_hat), abs(k-k_hat));
 
% Simulate using estimated parameters
t_end = t(end);
[t_hat, x1_hat, ~] = simulate_system(m_hat, k_hat, c_hat, t_end, dt);
 
% Compute position error
e_x = x1 - x1_hat;
 
% Plot: x(t) vs x_hat(t)
figure('Name', 'Exercise 2a — Least Squares Estimation', ...
       'NumberTitle', 'off', 'Position', [100, 100, 900, 600]);
 
subplot(2, 1, 1);
plot(t, x1, 'b', 'LineWidth', 1.5); hold on;
plot(t_hat, x1_hat, 'r--', 'LineWidth', 1.5);
xlabel('Time t [s]'); ylabel('Position [m]');
title('Position: x(t) vs $\hat{x}(t)$', 'Interpreter', 'latex');
legend('$x(t)$ — true', '$\hat{x}(t)$ — estimated', 'Interpreter', 'latex', 'Location', 'best');
grid on;
 
subplot(2, 1, 2);
plot(t, e_x, 'y', 'LineWidth', 1.5);
xlabel('Time t [s]'); ylabel('e_x(t) [m]');
title('Position Error: $e_x(t) = x(t) - \hat{x}(t)$', 'Interpreter', 'latex');
grid on;
 
sgtitle('Exercise 2a: Least Squares Estimation with T_s = 0.05 s', ...
        'FontSize', 13, 'FontWeight', 'bold');
 
%% =========================================================
%  Exercise 2b — Effect of Sampling Period Ts on Accuracy
%% =========================================================
fprintf('===== Exercise 2b =====\n');
 
Ts_values = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0];
 
error_m = zeros(size(Ts_values));
error_c = zeros(size(Ts_values));
error_k = zeros(size(Ts_values));
 
for i = 1 : length(Ts_values)
    Ts_i = Ts_values(i);
    step_i = round(Ts_i / dt);

    if step_i > length(t)
        continue; 
    end
 
    idx_i = 1 : step_i : length(t);
    x1_i = x1(idx_i);
    x2_i = x2(idx_i);
    u_i = u_out(idx_i);
 
    % Estimate acceleration
    x2_dot_i = gradient(x2_i, Ts_i);
 
    [m_i, c_i, k_i] = least_squares(x1_i, x2_i, x2_dot_i, u_i);
 
    error_m(i) = abs(m - m_i);
    error_c(i) = abs(c - c_i);
    error_k(i) = abs(k - k_i);
 
    fprintf('Ts = %.3f s → error_m=%.5f, error_c=%.5f, error_k=%.5f\n', ...
            Ts_i, error_m(i), error_c(i), error_k(i));
end
 
fprintf('\n');
 
% Plot: Parameter errors vs Ts
figure('Name', 'Exercise 2b — Effect of Ts', ...
       'NumberTitle', 'off', 'Position', [100, 100, 900, 600]);
 
subplot(3, 1, 1);
semilogx(Ts_values, error_m, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Sampling period T_s [s]'); ylabel('$|m - \hat{m}|$', 'Interpreter', 'latex');
title('Estimation error: mass m'); grid on;
 
subplot(3, 1, 2);
semilogx(Ts_values, error_c, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Sampling period T_s [s]'); ylabel('$|c - \hat{c}|$', 'Interpreter', 'latex');
title('Estimation error: damping c'); grid on;
 
subplot(3, 1, 3);
semilogx(Ts_values, error_k, 'go-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Sampling period T_s [s]'); ylabel('$|k - \hat{k}|$', 'Interpreter', 'latex');
title('Estimation error: spring constant k'); grid on;
 
sgtitle('Exercise 2b: Least Squares Estimation Error vs Sampling Period T_s', ...
        'FontSize', 13, 'FontWeight', 'bold');
 
%% =========================================================
%  Exercise 2c — Effect of Gaussian Noise on Estimation
%% =========================================================
fprintf('===== Exercise 2c =====\n');
 
Ts = 0.05;
step = round(Ts / dt);
idx = 1 : step : length(t);
t_samp  = t(idx);
x1_samp = x1(idx);
x2_samp = x2(idx);
u_samp  = u_out(idx);
N = length(x1_samp);
 
noise_levels = [0.001, 0.005, 0.01, 0.05];
 
% --- Clean estimate (no noise, reference) ---
x2_dot_clean = gradient(x2_samp, Ts);
 
[m_clean, c_clean, k_clean] = least_squares(x1_samp, x2_samp, x2_dot_clean, u_samp);
fprintf('No noise → m=%.4f, c=%.4f, k=%.4f\n', m_clean, c_clean, k_clean);
 
% --- Estimates with noise ---
rng(42);
 
err_m_n = zeros(size(noise_levels));
err_c_n = zeros(size(noise_levels));
err_k_n = zeros(size(noise_levels));
 
for i = 1 : length(noise_levels)
    sigma = noise_levels(i);
 
    x1_noisy = x1_samp + sigma * randn(N, 1);
    x2_noisy = x2_samp + sigma * randn(N, 1);
    u_noisy  = u_samp  + sigma * randn(N, 1);
 
    x2_dot_n = gradient(x2_noisy, Ts);
 
    [m_n, c_n, k_n] = least_squares(x1_noisy, x2_noisy, x2_dot_n, u_noisy);
 
    err_m_n(i) = abs(m - m_n);
    err_c_n(i) = abs(c - c_n);
    err_k_n(i) = abs(k - k_n);
 
    fprintf('sigma = %.3f → m=%.4f (err=%.4f), c=%.4f (err=%.4f), k=%.4f (err=%.4f)\n', ...
            sigma, m_n, err_m_n(i), c_n, err_c_n(i), k_n, err_k_n(i));
end
 
fprintf('\n');
 
% --- Plot: Parameter errors vs noise level ---
figure('Name', 'Theme 2c — Effect of Noise', ...
       'NumberTitle', 'off', 'Position', [100, 100, 900, 600]);
 
subplot(3, 1, 1);
semilogx(noise_levels, err_m_n, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Noise std \sigma'); ylabel('$|m - \hat{m}|$', 'Interpreter', 'latex');
title('Estimation error: mass m'); grid on;
 
subplot(3, 1, 2);
semilogx(noise_levels, err_c_n, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Noise std \sigma'); ylabel('$|c - \hat{c}|$', 'Interpreter', 'latex');
title('Estimation error: damping c'); grid on;
 
subplot(3, 1, 3);
semilogx(noise_levels, err_k_n, 'go-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Noise std \sigma'); ylabel('$|k - \hat{k}|$', 'Interpreter', 'latex');
title('Estimation error: spring constant k'); grid on;
 
sgtitle('Theme 2c: LS Estimation Error vs Noise Level \sigma', ...
        'FontSize', 13, 'FontWeight', 'bold');
 
% --- Plot: clean vs noisy signal comparison ---
figure('Name', 'Theme 2c — Clean vs Noisy', ...
       'NumberTitle', 'off', 'Position', [100, 100, 900, 400]);
 
sigma_demo    = noise_levels(end);
x1_noisy_demo = x1_samp + sigma_demo * randn(N, 1);
 
plot(t_samp, x1_samp,          'b',  'LineWidth', 1.5); hold on;
plot(t_samp, x1_noisy_demo, 'r.', 'MarkerSize', 4);
xlabel('Time t [s]'); ylabel('x(t) [m]');
title(sprintf('Clean vs Noisy position signal (\\sigma = %.3f)', sigma_demo));
legend('Clean x(t)', 'Noisy x(t)', 'Location', 'best');
grid on;
 
%% --- Save Results ---
if ~exist('output', 'dir')
    mkdir('output');
end
 
saveas(figure(1), 'output/tema2a_estimation.png');
saveas(figure(2), 'output/tema2b_sampling.png');
saveas(figure(3), 'output/tema2c_noise.png');
saveas(figure(4), 'output/tema2c_clean_vs_noisy.png');
 
save('output/tema2_data.mat', 'm_hat', 'c_hat', 'k_hat', ...
     'Ts_values', 'error_m', 'error_c', 'error_k', ...
     'noise_levels', 'err_m_n', 'err_c_n', 'err_k_n');
 
fprintf('Figures and data saved in /output folder.\n');
 