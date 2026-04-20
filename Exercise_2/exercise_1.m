%% Simulation and Modelling of Dynamic Systems
% Lab 02 — Exercise 1
% Gradient Method Parameter Estimation

clear; clc; close all;

%% System Parameters
J = 0.025;  % moment of inertia (known)
k = 0.30;   % damping coefficient k in [0.1, 0.5]
b = 1.00;   % input gain constant b in [0.5, 1.5]

fprintf('===== True System Parameters =====\n');
fprintf('J = %.4f kg·m² (known)\n', J);
fprintf('k = %.4f       (unknown to estimator)\n', k);
fprintf('b = %.4f       (unknown to estimator)\n\n', b);

%% Simulation Settings
dt = 1e-3;    % integration step [s]
t_end = 100;  % simulation duration [s]

%% Estimator Settings
fprintf('===== Gradient Method Estimation =====\n');
gamma_k = 50.0;  % adaptation gain for k
gamma_b = 150.0;  % adaptation gain for b
k0 = 0.10;       % initial estimate (lower bound of range)
b0 = 0.50;       % initial estimate (lower bound of range)

fprintf('Initial estimates: k0 = %.4f, b0 = %.4f\n', k0, b0);
fprintf('Adaptation gains:  gamma_k = %.4f, gamma_b = %.4f\n\n', gamma_k, gamma_b);

%% Exercise 1a — Gradient Method, No Disturbance d(t) = 0
fprintf('===== Exercise 1a — Gradient Method, No Disturbance =====\n');

d_none = @(t) 0;

fprintf('Simulating true system with d(t) = 0 ...\n');
[t, phi, phi_dot, phi_ddot] = simulate_vehicle(J, k, b, t_end, dt, d_none);
u = 0.25 * sin(0.5 * pi * t);
fprintf('Total time points: %d\n\n', length(t));

% Gradient estimation
[k_hat, b_hat] = gradient_estimator(t, phi_dot, phi_ddot, u, ...
                                     J, gamma_k, gamma_b, k0, b0);

% Simulate with final estimated parameters to obtain phi_hat
[t_hat, phi_hat, ~, ~] = simulate_vehicle(J, k_hat(end), b_hat(end), t_end, dt, d_none);

% Error signals
error_phi = phi - phi_hat;
error_k = k_hat - k;
error_b = b_hat - b;

fprintf('True parameters:      k = %.4f,      b = %.4f\n', k, b);
fprintf('Final estimates:  k_hat = %.4f,  b_hat = %.4f\n', k_hat(end), b_hat(end));
fprintf('Absolute errors:  |e_k| = %.5f, |e_b| = %.5f\n\n', ...
        abs(k - k_hat(end)), abs(b - b_hat(end)));

%% Persistent Excitation (PE) Metric Calculation
% Create the regressor vector
W = [-phi_dot, u];

% Calculate the excitation matrix R = (1/t_end) * Integral(W * W') dt
R_matrix = zeros(2, 2);
for i = 1 : length(t)-1
    dt_i = t(i+1) - t(i);
    w_i = W(i, :)';
    R_matrix = R_matrix + (w_i * w_i') * dt_i;
end
R_matrix = R_matrix / t_end;

% Calculate Eigenvalues
eigenvalues = eig(R_matrix);
lambda_min = min(eigenvalues);

fprintf('===== Persistent Excitation (PE) Analysis (d=0) =====\n');
fprintf('Eigenvalues of Excitation Matrix R: [%.6f,  %.6f]\n', eigenvalues(1), eigenvalues(2));
fprintf('Minimum Eigenvalue (PE Metric): lambda_min = %.6f\n\n', lambda_min);

%% Exercise 1a — Plots
% Plot 1a-1: True vs Estimated Angle
figure('Name', 'Exercise 1a - Angle Estimation', ...
       'NumberTitle', 'off', 'Position', [80, 100, 900, 600]);

subplot(2, 1, 1);
plot(t, phi, 'b', 'LineWidth', 1.5); hold on;
plot(t_hat, phi_hat, 'r--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\phi(t)$ [rad]', 'Interpreter', 'latex');
title('Orientation Angle: $\phi(t)$ vs $\hat{\phi}(t)$', 'Interpreter', 'latex');
legend('$\phi(t)$ - true', '$\hat{\phi}(t)$ - estimated', ...
       'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t, error_phi, 'm', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_\phi(t)$ [rad]', 'Interpreter', 'latex');
title('Angle Error: $e_\phi(t) = \phi(t) - \hat{\phi}(t)$', 'Interpreter', 'latex');
yline(0, 'k--', 'LineWidth', 0.8);
grid on;

sgtitle('Exercise 1a: Gradient Method - Angle Estimation ($d = 0$)', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Plot 1a-2: Parameter k Estimation and Error
figure('Name', 'Exercise 1a - Parameter k Estimation', ...
       'NumberTitle', 'off', 'Position', [100, 120, 900, 600]);
subplot(2, 1, 1);
plot(t, k_hat, 'b', 'LineWidth', 1.5); hold on;
yline(k, 'r--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\hat{k}(t)$', 'Interpreter', 'latex');
title('Parameter Estimate: $\hat{k}(t)$ vs True $k$', 'Interpreter', 'latex');
legend('$\hat{k}(t)$ - estimated', 'True $k$', 'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t, error_k, 'm', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_k(t)$', 'Interpreter', 'latex');
title('Parameter Error: $e_k(t) = \hat{k}(t) - k$', 'Interpreter', 'latex');
yline(0, 'k--', 'LineWidth', 0.8);
grid on;

sgtitle('Exercise 1a: Gradient Method - Parameter $k$ Estimation ($d = 0$)', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Plot 1a-3: Parameter b Estimation and Error
figure('Name', 'Exercise 1a - Parameter b Estimation', ...
       'NumberTitle', 'off', 'Position', [120, 140, 900, 600]);
subplot(2, 1, 1);
plot(t, b_hat, 'b', 'LineWidth', 1.5); hold on;
yline(b, 'r--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\hat{b}(t)$', 'Interpreter', 'latex');
title('Parameter Estimate: $\hat{b}(t)$ vs True $b$', 'Interpreter', 'latex');
legend('$\hat{b}(t)$ - estimated', 'True $b$', 'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t, error_b, 'm', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_b(t)$', 'Interpreter', 'latex');
title('Parameter Error: $e_b(t) = \hat{b}(t) - b$', 'Interpreter', 'latex');
yline(0, 'k--', 'LineWidth', 0.8);
grid on;

sgtitle('Exercise 1a: Gradient Method - Parameter $b$ Estimation ($d = 0$)', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

%% Exercise 1b — Gradient Method, With Disturbance d(t) = 0.02*sin(2t)
fprintf('===== Exercise 1b — Gradient Method, With Disturbance =====\n');

d_func = @(t) 0.02 * sin(2.0 * t);

fprintf('Simulating true system with d(t) = 0.02*sin(2.0*t) ...\n');
[t_d, phi_d, phi_dot_d, phi_ddot_d] = simulate_vehicle(J, k, b, t_end, dt, d_func);
u_d = 0.25 * sin(0.5 * pi * t_d);
fprintf('Total time points: %d\n\n', length(t_d));

% Gradient estimation on disturbed measurements
% Note: phi_ddot_d already contains d(t) effects — the estimator does NOT know d(t)
[k_hat_d, b_hat_d] = gradient_estimator(t_d, phi_dot_d, phi_ddot_d, u_d, ...
                                          J, gamma_k, gamma_b, k0, b0);

% Simulate with final estimated parameters (disturbance-free model for phi_hat)
[t_hat_d, phi_hat_d, ~, ~] = simulate_vehicle(J, k_hat_d(end), b_hat_d(end), t_end, dt, d_none);

% Error signals (with disturbance)
error_phi_d = phi_d - phi_hat_d;
error_k_d = k_hat_d - k;
error_b_d = b_hat_d - b;

fprintf('True parameters:           k = %.4f,      b = %.4f\n', k, b);
fprintf('Final estimates (d=0): k_hat = %.4f,  b_hat = %.4f\n', ...
        k_hat(end), b_hat(end));
fprintf('Final estimates (d≠0): k_hat = %.4f,  b_hat = %.4f\n', ...
        k_hat_d(end), b_hat_d(end));
fprintf('Absolute errors (d=0): |e_k| = %.5f, |e_b| = %.5f\n', ...
        abs(k - k_hat(end)),   abs(b - b_hat(end)));
fprintf('Absolute errors (d≠0): |e_k| = %.5f, |e_b| = %.5f\n\n', ...
        abs(k - k_hat_d(end)), abs(b - b_hat_d(end)));

%% Persistent Excitation (PE) Metric Calculation (With Disturbance)
% Create the regressor vector
W_d = [-phi_dot_d, u_d];

% Calculate the excitation matrix R_d = (1/t_end) * Integral(W_d * W_d') dt
R_matrix_d = zeros(2, 2);
for i = 1 : length(t_d)-1
    dt_i = t_d(i+1) - t_d(i);
    w_i = W_d(i, :)';
    R_matrix_d = R_matrix_d + (w_i * w_i') * dt_i;
end
R_matrix_d = R_matrix_d / t_end;

% Calculate Eigenvalues
eigenvalues_d = eig(R_matrix_d);
lambda_min_d = min(eigenvalues_d);

fprintf('===== Persistent Excitation (PE) Analysis (d≠0) =====\n');
fprintf('Eigenvalues of Excitation Matrix R: [%.6f,  %.6f]\n', eigenvalues_d(1), eigenvalues_d(2));
fprintf('Minimum Eigenvalue (PE Metric): lambda_min = %.6f\n\n', lambda_min_d);

%% Exercise 1b - Plots
% Plot 1b-1: True vs Estimated Angle (With Disturbance)
figure('Name', 'Exercise 1b - Angle Estimation', ...
       'NumberTitle', 'off', 'Position', [80, 100, 900, 600]);
subplot(2, 1, 1);
plot(t_d, phi_d, 'b', 'LineWidth', 1.5); hold on;
plot(t_hat_d, phi_hat_d, 'r--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\phi(t)$ [rad]', 'Interpreter', 'latex');
title('Orientation Angle: $\phi(t)$ vs $\hat{\phi}(t)$ (with $d \neq 0$)', 'Interpreter', 'latex');
legend('$\phi(t)$ - true', '$\hat{\phi}(t)$ - estimated', ...
       'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t_d, error_phi_d, 'm', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_\phi(t)$ [rad]', 'Interpreter', 'latex');
title('Angle Error: $e_\phi(t) = \phi(t) - \hat{\phi}(t)$ (with $d\neq 0$)', 'Interpreter', 'latex');
yline(0, 'k--', 'LineWidth', 0.8);
grid on;

sgtitle('Exercise 1b: Gradient Method - Angle Estimation ($d(t) \neq 0$)', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Plot 1b-2: Parameter k Estimation and Error (With Disturbance)
figure('Name', 'Exercise 1b - Parameter k Estimation', ...
       'NumberTitle', 'off', 'Position', [100, 120, 900, 600]);
subplot(2, 1, 1);
plot(t_d, k_hat_d, 'b', 'LineWidth', 1.5); hold on;
yline(k, 'r--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\hat{k}(t)$', 'Interpreter', 'latex');
title('Parameter Estimate: $\hat{k}(t)$ vs True $k$ (with $d \neq 0$)', 'Interpreter', 'latex');
legend('$\hat{k}(t)$ - estimated', 'True $k$', 'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t_d, error_k_d, 'm', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_k(t)$', 'Interpreter', 'latex');
title('Parameter Error: $e_k(t) = \hat{k}(t) - k$', 'Interpreter', 'latex');
yline(0, 'k--', 'LineWidth', 0.8);
grid on;

sgtitle('Exercise 1b: Gradient Method - Parameter $k$ Estimation ($d(t) \neq 0$)', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Plot 1b-3: Parameter b Estimation and Error (With Disturbance)
figure('Name', 'Exercise 1b - Parameter b Estimation', ...
       'NumberTitle', 'off', 'Position', [120, 140, 900, 600]);
subplot(2, 1, 1);
plot(t_d, b_hat_d, 'b', 'LineWidth', 1.5); hold on;
yline(b, 'r--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\hat{b}(t)$', 'Interpreter', 'latex');
title('Parameter Estimate: $\hat{b}(t)$ vs True $b$ (with $d \neq 0$)', 'Interpreter', 'latex');
legend('$\hat{b}(t)$ - estimated', 'True $b$', 'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t_d, error_b_d, 'm', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_b(t)$', 'Interpreter', 'latex');
title('Parameter Error: $e_b(t) = \hat{b}(t) - b$', 'Interpreter', 'latex');
yline(0, 'k--', 'LineWidth', 0.8);
grid on;

sgtitle('Exercise 1b: Gradient Method - Parameter $b$ Estimation ($d(t) \neq 0$)', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

%% Plot 1b-4: Comparison (d=0 vs d!=0)
figure('Name', 'Exercise 1 - Disturbance Comparison', ...
       'NumberTitle', 'off', 'Position', [140, 20, 800, 800]);

% Subplot 4.1: phi_hat comparison
subplot(4, 1, 1);
plot(t, phi_hat, 'b', 'LineWidth', 1.5); hold on;
plot(t_d, phi_hat_d, 'r--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\hat{\phi}(t)$ [rad]', 'Interpreter', 'latex');
title('Estimated Angle $\hat{\phi}(t)$ Comparison', 'Interpreter', 'latex');
legend('Clean ($d=0$)', 'Disturbed ($d \neq 0$)', 'Interpreter', 'latex', 'Location', 'best');
grid on;

% Subplot 4.2: e_phi comparison
subplot(4, 1, 2);
plot(t, error_phi, 'b', 'LineWidth', 1.5); hold on;
plot(t_d, error_phi_d, 'r', 'LineWidth', 1.5);
yline(0, 'k--', 'LineWidth', 0.8);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_\phi(t)$ [rad]', 'Interpreter', 'latex');
title('Angle Error $e_\phi(t)$ Comparison', 'Interpreter', 'latex');
legend('Clean ($d=0$)', 'Disturbed ($d \neq 0$)', 'Interpreter', 'latex', 'Location', 'best');
grid on;

% Subplot 4.3: e_k comparison
subplot(4, 1, 3);
plot(t, error_k, 'b', 'LineWidth', 1.5); hold on;
plot(t_d, error_k_d, 'r', 'LineWidth', 1.5);
yline(0, 'k--', 'LineWidth', 0.8);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_k(t)$', 'Interpreter', 'latex');
title('Parameter Error $e_k(t)$ Comparison', 'Interpreter', 'latex');
legend('Clean ($d=0$)', 'Disturbed ($d \neq 0$)', 'Interpreter', 'latex', 'Location', 'best');
grid on;

% Subplot 4.4: e_b comparison
subplot(4, 1, 4);
plot(t, error_b, 'b', 'LineWidth', 1.5); hold on;
plot(t_d, error_b_d, 'r', 'LineWidth', 1.5);
yline(0, 'k--', 'LineWidth', 0.8);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_b(t)$', 'Interpreter', 'latex');
title('Parameter Error $e_b(t)$ Comparison', 'Interpreter', 'latex');
legend('Clean ($d=0$)', 'Disturbed ($d \neq 0$)', 'Interpreter', 'latex', 'Location', 'best');
grid on;

sgtitle('Exercise 1: System Performance Comparison ($d=0$ vs $d \neq 0$)', ...
        'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');

%% Summary
fprintf('===== Summary =====\n');
fprintf('Case d=0: k_hat=%.5f (|e_k|=%.5f), b_hat=%.5f (|e_b|=%.5f)\n', ...
        k_hat(end), abs(k - k_hat(end)), b_hat(end), abs(b - b_hat(end)));
fprintf('Case d≠0: k_hat=%.5f (|e_k|=%.5f), b_hat=%.5f (|e_b|=%.5f)\n', ...
        k_hat_d(end), abs(k - k_hat_d(end)), b_hat_d(end), abs(b - b_hat_d(end)));