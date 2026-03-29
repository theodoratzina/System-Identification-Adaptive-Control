%% Simulation and Modelling of Dynamic Systems
% Lab 01 - Exercise 2

clear; clc; close all;

%% Load data from Exercise 1
% Loads: t, x1, x2, u_out, m, k, c, dt
load('simulation_data.mat');
 
fprintf('===== Data loaded from simulation_data.mat =====\n');
fprintf('True parameters: m=%.4f, k=%.4f, c=%.4f\n\n', m, k, c);

%%  Exercise 2a — Least Squares Estimation (Ts = 0.05 s)
fprintf('===== Exercise 2a =====\n');
Ts = 0.05;  % sampling period [s]
 
% Sample the simulation data at every Ts seconds
step = round(Ts / dt);
idx  = 1 : step : length(t);
t_samp = t(idx);
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
xlabel('Time $t$ [s]', 'Interpreter', 'latex'); 
ylabel('Position [m]', 'Interpreter', 'latex');
title('Position: $x(t)$ vs $\hat{x}(t)$', 'Interpreter', 'latex');
legend('$x(t)$ — true', '$\hat{x}(t)$ — estimated', 'Interpreter', 'latex', 'Location', 'best');
grid on;
 
subplot(2, 1, 2);
plot(t, e_x, 'y', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex'); 
ylabel('$e_x(t)$ [m]', 'Interpreter', 'latex');
title('Position Error: $e_x(t) = x(t) - \hat{x}(t)$', 'Interpreter', 'latex');
grid on;
 
sgtitle('Exercise 2a: Least Squares Estimation with $T_s = 0.05$ s', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');
 
%%  Exercise 2b — Effect of Sampling Period Ts on Accuracy
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
xlabel('Sampling period $T_s$ [s]', 'Interpreter', 'latex'); 
ylabel('$|m - \hat{m}|$', 'Interpreter', 'latex');
title('Estimation error: mass $m$', 'Interpreter', 'latex');
grid on;
 
subplot(3, 1, 2);
semilogx(Ts_values, error_c, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Sampling period $T_s$ [s]', 'Interpreter', 'latex'); 
ylabel('$|c - \hat{c}|$', 'Interpreter', 'latex');
title('Estimation error: damping $c$', 'Interpreter', 'latex');
grid on;
 
subplot(3, 1, 3);
semilogx(Ts_values, error_k, 'go-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Sampling period $T_s$ [s]', 'Interpreter', 'latex'); 
ylabel('$|k - \hat{k}|$', 'Interpreter', 'latex');
title('Estimation error: spring constant $k$', 'Interpreter', 'latex'); 
grid on;
 
sgtitle('Exercise 2b: Least Squares Estimation Error vs Sampling Period $T_s$', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');
 
%%  Exercise 2c — Effect of Gaussian Noise on Estimation
% Reuses t_samp, x1_samp, x2_samp, u_samp, x2_dot from 2a (same Ts = 0.05 s)
% Input u(t) is known analytically — no noise added
fprintf('===== Exercise 2c =====\n');
 
% Noise levels as percentage of signal standard deviation
noise_values = [0.001, 0.005, 0.01, 0.05, 0.1];  % 0.1% to 10%
sigma_x = noise_values * std(x1_samp);
sigma_v = noise_values * std(x2_samp);
 
error_m_noisy = zeros(size(noise_values));
error_c_noisy = zeros(size(noise_values));
error_k_noisy = zeros(size(noise_values));
 
rng(42);
fprintf('True values: m=%.4f, c=%.4f, k=%.4f\n', m, c, k);
fprintf('No noise   : m=%.4f, c=%.4f, k=%.4f\n\n', m_hat, c_hat, k_hat);
 
for i = 1 : length(noise_values)
 
    % Add white Gaussian noise to position and velocity
    x1_noisy = x1_samp + sigma_x(i) * randn(size(x1_samp));
    x2_noisy = x2_samp + sigma_v(i) * randn(size(x2_samp));
 
    % Estimate acceleration from noisy velocity
    x2_dot_noisy = gradient(x2_noisy, Ts);
 
    % Least squares on noisy data
    [m_noisy, c_noisy, k_noisy] = least_squares(x1_noisy, x2_noisy, x2_dot_noisy, u_samp);
 
    error_m_noisy(i) = abs(m - m_noisy);
    error_c_noisy(i) = abs(c - c_noisy);
    error_k_noisy(i) = abs(k - k_noisy);
 
    fprintf('Noise=%.1f%%  ->  m=%.4f (e=%.5f), c=%.4f (e=%.5f), k=%.4f (e=%.5f)\n', ...
            noise_values(i)*100, m_noisy, error_m_noisy(i), c_noisy, error_c_noisy(i), k_noisy, error_k_noisy(i));
end
 
fprintf('\n');
 
% Plot: parameter estimation errors vs noise level
figure('Name', 'Exercise 2c — Effect of Noise', ...
       'NumberTitle', 'off', 'Position', [100, 100, 900, 600]);
 
subplot(3, 1, 1);
semilogx(noise_values*100, error_m_noisy, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Noise level [\% of signal std]', 'Interpreter', 'latex');
ylabel('$|m - \hat{m}|$', 'Interpreter', 'latex');
title('Estimation error: mass $m$', 'Interpreter', 'latex');
grid on;
 
subplot(3, 1, 2);
semilogx(noise_values*100, error_c_noisy, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Noise level [\% of signal std]', 'Interpreter', 'latex');
ylabel('$|c - \hat{c}|$', 'Interpreter', 'latex');
title('Estimation error: damping $c$', 'Interpreter', 'latex');
grid on;
 
subplot(3, 1, 3);
semilogx(noise_values*100, error_k_noisy, 'go-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Noise level [\% of signal std]', 'Interpreter', 'latex');
ylabel('$|k - \hat{k}|$', 'Interpreter', 'latex');
title('Estimation error: spring constant $k$', 'Interpreter', 'latex');
grid on;
 
sgtitle('Exercise 2c: Least Squares Estimation Error vs Noise Level', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');
 
% Plot: clean vs noisy signal (highest noise level)
figure('Name', 'Exercise 2c — Clean vs Noisy', ...
       'NumberTitle', 'off', 'Position', [100, 100, 900, 400]);
 
x1_noisy_demo = x1_samp + sigma_x(end) * randn(size(x1_samp));
plot(t_samp, x1_samp, 'b',  'LineWidth', 1.5); hold on;
plot(t_samp, x1_noisy_demo, 'r.', 'MarkerSize', 8);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$x(t)$ [m]', 'Interpreter', 'latex');
title(sprintf('Clean vs Noisy position (noise = %.0f\\%% of std)', ...
    noise_values(end)*100), 'Interpreter', 'latex');
legend('Clean $x(t)$', 'Noisy $x(t)$', 'Interpreter', 'latex', 'Location', 'best');
grid on;