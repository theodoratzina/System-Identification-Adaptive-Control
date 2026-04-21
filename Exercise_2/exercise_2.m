%% Simulation and Modelling of Dynamic Systems
% Lab 02 — Exercise 2
% Lyapunov-Based Parameter Estimation — Simple Pendulum

clear; clc; close all;

%% System Parameters
g = 9.81;  % gravitational acceleration [m/s²] (known)
l = 1.80;  % pendulum length [m] l in [1.0, 2.5]
c = 0.50;  % damping coefficient c in [0.2, 0.8]
alpha = g / l;

fprintf('===== True System Parameters =====\n');
fprintf('g = %.4f m/s² (known)\n', g);
fprintf('l = %.4f m    (unknown to estimator)\n', l);
fprintf('c = %.4f      (unknown to estimator)\n', c);
fprintf('alpha = g/l = %.4f rad/s²\n\n', alpha);

%% Simulation Settings
dt = 1e-3;   % integration step [s]
t_end = 100;  % simulation duration [s]

%% Estimator Settings
gamma_alpha = 70.0;
gamma_c = 15.0;
alpha0 = g / 1.0;  % assume l = 1.0 m
c0 = 0.20;         % initial estimate (lower bound of range)

fprintf('Initial estimates: l0 = 1.00 m  (alpha0 = %.4f),  c0 = %.4f\n', alpha0, c0);
fprintf('Adaptation gains:  gamma_alpha = %.4f, gamma_c = %.4f\n\n', gamma_alpha, gamma_c);

%% Simulate True Pendulum
fprintf('Simulating true pendulum ...\n');
[t, theta, theta_dot] = simulate_pendulum(g, l, c, t_end, dt);
u = 0.5 * sin(t);
fprintf('Total time points: %d\n\n', length(t));

%% Exercise 2a — Lyapunov Method, No Noise
fprintf('===== Exercise 2a — Lyapunov Method, No Noise =====\n');

% Lyapunov estimation
[alpha_hat, c_hat, theta_hat, theta_dot_hat] = lyapunov_estimator( ...
    t, theta, theta_dot, u, gamma_alpha, gamma_c, alpha0, c0);

% Recover estimated pendulum length from estimated alpha
l_hat = g ./ alpha_hat;

% Error signals
error_theta = theta - theta_hat;
error_l = l_hat - l;
error_c = c_hat - c;

% Lyapunov driving signal
s = theta_dot - theta_dot_hat;

fprintf('True values:         l = %.4f m,     c = %.4f\n', l, c);
fprintf('Final estimates: l_hat = %.4f m, c_hat = %.4f\n', l_hat(end), c_hat(end));
fprintf('Absolute errors: |e_l| = %.5f,  |e_c| = %.5f\n\n', ...
        abs(l - l_hat(end)), abs(c - c_hat(end)));

%% Exercise 2a — Plots
% Plot 2a-1: True vs Estimated Angle
figure('Name', 'Exercise 2a - Angle Estimation', ...
       'NumberTitle', 'off', 'Position', [80, 100, 900, 600]);

subplot(2, 1, 1);
plot(t, theta, 'b', 'LineWidth', 1.5); hold on;
plot(t, theta_hat, 'r--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\theta(t)$ [rad]', 'Interpreter', 'latex');
title('Pendulum Angle: $\theta(t)$ vs $\hat{\theta}(t)$', 'Interpreter', 'latex');
legend('$\theta(t)$ - true', '$\hat{\theta}(t)$ - estimated', ...
       'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t, error_theta, 'm', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_\theta(t)$ [rad]', 'Interpreter', 'latex');
title('Angle Tracking Error: $e_\theta(t) = \theta(t) - \hat{\theta}(t)$', 'Interpreter', 'latex');
yline(0, 'k--', 'LineWidth', 0.8);
grid on;

sgtitle('Exercise 2a: Lyapunov Method - Angle Estimation (no noise)', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Plot 2a-2: Parameter l Estimation and Error
figure('Name', 'Exercise 2a - Parameter l Estimation', ...
       'NumberTitle', 'off', 'Position', [100, 120, 900, 600]);

subplot(2, 1, 1);
plot(t, l_hat, 'b', 'LineWidth', 1.5); hold on;
yline(l, 'r--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\hat{l}(t)$ [m]', 'Interpreter', 'latex');
title('Parameter Estimate: $\hat{l}(t) = g / \hat{\alpha}(t)$ vs True $l$', 'Interpreter', 'latex');
legend('$\hat{l}(t)$ - estimated', 'True $l$', 'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t, error_l, 'm', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_l(t)$ [m]', 'Interpreter', 'latex');
title('Parameter Error: $e_l(t) = \hat{l}(t) - l$', 'Interpreter', 'latex');
yline(0, 'k--', 'LineWidth', 0.8);
grid on;

sgtitle('Exercise 2a: Lyapunov Method - Parameter $l$ Estimation (no noise)', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Plot 2a-3: Parameter c Estimation and Error
figure('Name', 'Exercise 2a - Parameter c Estimation', ...
       'NumberTitle', 'off', 'Position', [120, 140, 900, 600]);

subplot(2, 1, 1);
plot(t, c_hat, 'b', 'LineWidth', 1.5); hold on;
yline(c, 'r--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\hat{c}(t)$', 'Interpreter', 'latex');
title('Parameter Estimate: $\hat{c}(t)$ vs True $c$', 'Interpreter', 'latex');
legend('$\hat{c}(t)$ - estimated', 'True $c$', 'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t, error_c, 'm', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_c(t)$', 'Interpreter', 'latex');
title('Parameter Error: $e_c(t) = \hat{c}(t) - c$', 'Interpreter', 'latex');
yline(0, 'k--', 'LineWidth', 0.8);
grid on;

sgtitle('Exercise 2a: Lyapunov Method - Parameter $c$ Estimation (no noise)', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Plot 2a-4: Lyapunov Driving Signal s(t)
figure('Name', 'Exercise 2a - Lyapunov Driving Signal', ...
       'NumberTitle', 'off', 'Position', [140, 160, 900, 400]);

plot(t, s, 'b', 'LineWidth', 1.5); hold on;
yline(0, 'r--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$s(t)$ [rad/s]', 'Interpreter', 'latex');
title({'Exercise 2a: Lyapunov Method - Velocity Tracking Error (no noise)'; ...
       'Lyapunov Driving Signal: $s(t) = \dot{\theta}(t) - \hat{\dot{\theta}}(t)$'}, ...
      'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
legend('$s(t)$ - driving signal', 'Zero reference', ...
       'Interpreter', 'latex', 'Location', 'best');
grid on;

%% Exercise 2b — Lyapunov Method, with Noise n(t) = n0 * sin(20*pi*t)
fprintf('===== Exercise 2b — Lyapunov Method, Noise Effect Study =====\n');
fprintf('Noise signal: n(t) = n0 * sin(20*pi*t)\n');
fprintf('Noise applied to: theta(t), theta_dot(t), u(t)\n\n');

n0_values = [0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2];

error_l = zeros(size(n0_values));
error_c = zeros(size(n0_values));

fprintf('%-12s %-12s %-12s %-12s %-12s\n', 'n0', 'l_hat', '|e_l|', 'c_hat', '|e_c|');
fprintf('%s\n', repmat('-', 1, 62));

for i = 1 : length(n0_values)
    n0 = n0_values(i);
    n = n0 * sin(20 * pi * t);

    theta_n = theta + n;
    theta_dot_n = theta_dot + n;
    u_n = u + n;

    [alpha_hat_i, c_hat_i, ~, ~] = lyapunov_estimator( ...
        t, theta_n, theta_dot_n, u_n, gamma_alpha, gamma_c, alpha0, c0);

    l_hat_final = g / alpha_hat_i(end);
    c_hat_final = c_hat_i(end);

    error_l(i) = abs(l - l_hat_final);
    error_c(i) = abs(c - c_hat_final);

    fprintf('%-12.4f %-12.5f %-12.5f %-12.5f %-12.5f\n', ...
            n0, l_hat_final, error_l(i), c_hat_final, error_c(i));
end

fprintf('\n');

%% Exercise 2b — Plots
% Define noise levels for the time-domain study
n0_low  = 0.01;
n0_high = 0.05;

% Run estimator for low noise
n_low = n0_low * sin(20 * pi * t);
[alpha_low, c_low, theta_hat_low, ~] = lyapunov_estimator( ...
    t, theta + n_low, theta_dot + n_low, u + n_low, gamma_alpha, gamma_c, alpha0, c0);
l_low = g ./ alpha_low;

% Run estimator for high noise
n_high = n0_high * sin(20 * pi * t);
[alpha_high, c_high, theta_hat_high, ~] = lyapunov_estimator( ...
    t, theta + n_high, theta_dot + n_high, u + n_high, gamma_alpha, gamma_c, alpha0, c0);
l_high = g ./ alpha_high;

% Plot 2b-1: Angle Tracking Error (Clean vs Low vs High)
figure('Name', 'Exercise 2b - Angle Tracking Comparison', ...
       'NumberTitle', 'off', 'Position', [100, 100, 900, 600]);

subplot(2, 1, 1);
plot(t, theta, 'g', 'LineWidth', 1.5); hold on;
plot(t, theta_hat_low, 'b--', 'LineWidth', 1.5);
plot(t, theta_hat_high, 'r:', 'LineWidth', 1.8);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\hat{\theta}(t)$ [rad]', 'Interpreter', 'latex');
title('Estimated Angle $\hat{\theta}(t)$ under Different Noise Levels', 'Interpreter', 'latex');
legend('True $\theta(t)$', sprintf('Low Noise ($\\eta_0 = %.2f$)', n0_low), ...
       sprintf('High Noise ($\\eta_0 = %.2f$)', n0_high), 'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t, error_theta, 'g', 'LineWidth', 1.5); hold on;
plot(t, theta - theta_hat_low, 'b--', 'LineWidth', 1.5);
plot(t, theta - theta_hat_high, 'r:', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_\theta(t)$ [rad]', 'Interpreter', 'latex');
title('Angle Tracking Error $e_\theta(t)$ under Different Noise Levels', 'Interpreter', 'latex');
legend('Clean ($\eta_0 = 0$)', sprintf('Low Noise ($\\eta_0 = %.2f$)', n0_low), ...
       sprintf('High Noise ($\\eta_0 = %.2f$)', n0_high), 'Interpreter', 'latex', 'Location', 'best');
grid on;

sgtitle('Exercise 2b: Angle Tracking Robustness vs Measurement Noise', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Plot 2b-2: Parameter Convergence (Clean vs Low vs High)
figure('Name', 'Exercise 2b - Parameter Convergence Comparison', ...
       'NumberTitle', 'off', 'Position', [120, 120, 900, 600]);

subplot(2, 1, 1);
plot(t, l_hat, 'm', 'LineWidth', 1.5); hold on;
plot(t, l_low, 'b--', 'LineWidth', 1.5);
plot(t, l_high, 'r:', 'LineWidth', 1.8);
yline(l, 'g-', 'LineWidth', 1.2);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\hat{l}(t)$ [m]', 'Interpreter', 'latex');
title('Pendulum Length Estimate $\hat{l}(t)$ for Different Noise Levels', 'Interpreter', 'latex');
legend('Clean ($\eta_0 = 0$)', sprintf('Low Noise ($\\eta_0 = %.2f$)', n0_low), ...
       sprintf('High Noise ($\\eta_0 = %.2f$)', n0_high), 'True $l$', ...
       'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t, c_hat, 'm', 'LineWidth', 1.5); hold on;
plot(t, c_low, 'b--', 'LineWidth', 1.5);
plot(t, c_high, 'r:', 'LineWidth', 1.8);
yline(c, 'g-', 'LineWidth', 1.2);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\hat{c}(t)$', 'Interpreter', 'latex');
title('Damping Estimate $\hat{c}(t)$ for Different Noise Levels', 'Interpreter', 'latex');
legend('Clean ($\eta_0 = 0$)', sprintf('Low Noise ($\\eta_0 = %.2f$)', n0_low), ...
       sprintf('High Noise ($\\eta_0 = %.2f$)', n0_high), 'True $c$', ...
       'Interpreter', 'latex', 'Location', 'best');
grid on;

sgtitle('Exercise 2b: Parameter Estimation Robustness vs Measurement Noise', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Plot 2b-3: Parameter Errors vs Noise Amplitude
figure('Name', 'Exercise 2b - Errors vs Noise Amplitude', ...
       'NumberTitle', 'off', 'Position', [140, 140, 900, 600]);

subplot(2, 1, 1);
semilogx(n0_values(2:end), error_l(2:end), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 7); hold on;
yline(error_l(1), 'r--', 'LineWidth', 1.5);
xlabel('Noise amplitude $\eta_0$', 'Interpreter', 'latex');
ylabel('$|l - \hat{l}|$ [m]', 'Interpreter', 'latex');
title('Final Estimation Error: Pendulum Length $l$', 'Interpreter', 'latex');
legend('With noise', 'Noise-free ($\eta_0 = 0$)', 'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
semilogx(n0_values(2:end), error_c(2:end), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 7); hold on;
yline(error_c(1), 'b--', 'LineWidth', 1.5);
xlabel('Noise amplitude $\eta_0$', 'Interpreter', 'latex');
ylabel('$|c - \hat{c}|$', 'Interpreter', 'latex');
title('Final Estimation Error: Damping Coefficient $c$', 'Interpreter', 'latex');
legend('With noise', 'Noise-free ($\eta_0 = 0$)', 'Interpreter', 'latex', 'Location', 'best');
grid on;

sgtitle('Exercise 2b: Final Parameter Estimation Error vs Noise Amplitude', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

%% ===== Summary =====
fprintf('===== Summary =====\n');
fprintf('%-12s %-12s %-12s\n', 'eta0', '|el|', '|ec|');
fprintf('%s\n', repmat('-', 1, 38));
for i = 1 : length(n0_values)
    fprintf('%-12.4f %-12.5f %-12.5f\n', n0_values(i), error_l(i), error_c(i));
end
