%% Simulation and Modelling of Dynamic Systems
% Lab 02 — Exercise 1
% Gradient Method Parameter Estimation

clear; clc; close all;

%% System Parameters
J = 0.025;  % moment of inertia (known)
k = 0.30;   % damping coefficient k in [0.1, 0.5]
b = 1.00;   % input gain constant b in [0.5, 1.5]

fprintf('===== True System Parameters =====\n');
fprintf('J = %.4f kg·m²  (known)\n', J);
fprintf('k = %.4f        (unknown to estimator)\n', k);
fprintf('b = %.4f        (unknown to estimator)\n\n', b);

%% Simulation Settings
dt = 1e-3;   % integration step [s]
t_end = 30;  % simulation duration [s]

%% Estimator Settings
fprintf('===== Gradient Method Estimation =====\n');
gamma_k = 0.50;  % adaptation gain for k
gamma_b = 0.50;  % adaptation gain for b
k0 = 0.10;       % initial estimate (lower bound of range)
b0 = 0.50;       % initial estimate (lower bound of range)

fprintf('Initial estimates: k0 = %.4f, b0 = %.4f\n\n', k0, b0);
fprintf('Adaptation gains:  gamma_k = %.4f, gamma_b = %.4f\n\n', gamma_k, gamma_b);

%% Exercise 1a — Gradient Method, No Disturbance d(t) = 0
fprintf('===== Exercise 1a — Gradient Method, No Disturbance =====\n');

d_none = @(t) 0;

fprintf('Simulating true system with d(t) = 0 ...\n');
[t, phi, phi_dot, phi_ddot] = simulate_vehicle(J, k, b, t_end, dt, d_none);
u_vec = 0.25 * sin(0.5 * pi * t);
fprintf('Total time points: %d\n\n', length(t));

% Gradient estimation
[k_hat, b_hat] = gradient_estimator(t, phi_dot, phi_ddot, u_vec, ...
                                     J, gamma_k, gamma_b, k0, b0);

% Simulate with final estimated parameters to obtain phi_hat
[t_hat, phi_hat, ~, ~] = simulate_vehicle(J, k_hat(end), b_hat(end), t_end, dt, d_none);

% Error signals
e_phi = phi - phi_hat;
e_k = k_hat - k;
e_b = b_hat - b;

fprintf('True parameters:      k = %.4f,     b = %.4f\n', k, b);
fprintf('Final estimates:  k_hat = %.4f, b_hat = %.4f\n', k_hat(end), b_hat(end));
fprintf('Absolute errors:  |e_k| = %.5f, |e_b| = %.5f\n\n', ...
        abs(k - k_hat(end)), abs(b - b_hat(end)));

%% Exercise 1a — Plots
% Plot 1: True vs Estimated Angle
figure('Name', 'Exercise 1a — Angle Estimation', ...
       'NumberTitle', 'off', 'Position', [80, 100, 900, 600]);

subplot(2, 1, 1);
plot(t, phi, 'b', 'LineWidth', 1.5); hold on;
plot(t_hat, phi_hat, 'r--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\phi(t)$ [rad]', 'Interpreter', 'latex');
title('Orientation Angle: $\phi(t)$ vs $\hat{\phi}(t)$', 'Interpreter', 'latex');
legend('$\phi(t)$ — true', '$\hat{\phi}(t)$ — estimated', ...
       'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t, e_phi, 'm', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_\phi(t)$ [rad]', 'Interpreter', 'latex');
title('Angle Error: $e_\phi(t) = \phi(t) - \hat{\phi}(t)$', 'Interpreter', 'latex');
yline(0, 'k--', 'LineWidth', 0.8);
grid on;

sgtitle('Exercise 1a: Gradient Method — Angle Estimation ($d = 0$)', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Plot 2: Parameter Estimate Convergence
figure('Name', 'Exercise 1a — Parameter Convergence', ...
       'NumberTitle', 'off', 'Position', [80, 100, 900, 600]);

subplot(2, 1, 1);
plot(t, k_hat, 'b', 'LineWidth', 1.5); hold on;
yline(k, 'r--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\hat{k}(t)$', 'Interpreter', 'latex');
title('Parameter Estimate: $\hat{k}(t)$', 'Interpreter', 'latex');
legend('$\hat{k}(t)$', 'True $k$', 'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t, b_hat, 'r', 'LineWidth', 1.5); hold on;
yline(b, 'b--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\hat{b}(t)$', 'Interpreter', 'latex');
title('Parameter Estimate: $\hat{b}(t)$', 'Interpreter', 'latex');
legend('$\hat{b}(t)$', 'True $b$', 'Interpreter', 'latex', 'Location', 'best');
grid on;

sgtitle('Exercise 1a: Gradient Method — Parameter Convergence ($d = 0$)', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Plot 3: Parameter Estimation Errors
figure('Name', 'Exercise 1a — Parameter Errors', ...
       'NumberTitle', 'off', 'Position', [80, 100, 900, 600]);

subplot(2, 1, 1);
plot(t, e_k, 'b', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_k(t) = \hat{k}(t) - k$', 'Interpreter', 'latex');
title('Parameter Error: $e_k(t)$', 'Interpreter', 'latex');
yline(0, 'k--', 'LineWidth', 0.8);
grid on;

subplot(2, 1, 2);
plot(t, e_b, 'r', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_b(t) = \hat{b}(t) - b$', 'Interpreter', 'latex');
title('Parameter Error: $e_b(t)$', 'Interpreter', 'latex');
yline(0, 'k--', 'LineWidth', 0.8);
grid on;

sgtitle('Exercise 1a: Gradient Method — Parameter Estimation Errors ($d = 0$)', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

%% Exercise 1b — Gradient Method, With Disturbance d(t) = 0.02*sin(2t)
fprintf('===== Exercise 1b — Gradient Method, With Disturbance =====\n');

d_func = @(t) 0.02 * sin(2.0 * t);

fprintf('Simulating true system with d(t) = 0.02*sin(2.0*t) ...\n');
[t_d, phi_d, phi_dot_d, phi_ddot_d] = simulate_vehicle(J, k, b, t_end, dt, d_func);
u_vec_d = 0.25 * sin(0.5 * pi * t_d);
fprintf('Total time points: %d\n\n', length(t_d));

% Gradient estimation on disturbed measurements
% Note: phi_ddot_d already contains d(t) effects — the estimator does NOT know d(t)
[k_hat_d, b_hat_d] = gradient_estimator(t_d, phi_dot_d, phi_ddot_d, u_vec_d, ...
                                          J, gamma_k, gamma_b, k0, b0);

% Simulate with final estimated parameters (disturbance-free model for phi_hat)
[t_hat_d, phi_hat_d, ~, ~] = simulate_vehicle(J, k_hat_d(end), b_hat_d(end), t_end, dt, d_none);

% Error signals (with disturbance)
e_phi_d = phi_d - phi_hat_d;
e_k_d   = k_hat_d - k;
e_b_d   = b_hat_d - b;

fprintf('True parameters:           k = %.4f,     b = %.4f\n', k, b);
fprintf('Final estimates (d=0): k_hat = %.4f, b_hat = %.4f\n', ...
        k_hat(end), b_hat(end));
fprintf('Final estimates (d≠0): k_hat = %.4f, b_hat = %.4f\n', ...
        k_hat_d(end), b_hat_d(end));
fprintf('Absolute errors (d=0): |e_k| = %.5f, |e_b| = %.5f\n', ...
        abs(k - k_hat(end)),   abs(b - b_hat(end)));
fprintf('Absolute errors (d≠0): |e_k| = %.5f, |e_b| = %.5f\n\n', ...
        abs(k - k_hat_d(end)), abs(b - b_hat_d(end)));

%% Exercise 1b — Plots
% Plot 1: Estimated vs True Angle (with disturbance)
figure('Name', 'Exercise 1b — Angle Estimation with Disturbance', ...
       'NumberTitle', 'off', 'Position', [80, 100, 900, 600]);

subplot(2, 1, 1);
plot(t_d, phi_d, 'b', 'LineWidth', 1.5); hold on;
plot(t_hat_d, phi_hat_d, 'r--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\phi(t)$ [rad]', 'Interpreter', 'latex');
title('Orientation Angle: $\phi(t)$ vs $\hat{\phi}(t)$ (with $d \neq 0$)', 'Interpreter', 'latex');
legend('$\phi(t)$ — true', '$\hat{\phi}(t)$ — estimated', ...
       'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t_d, e_phi_d, 'm', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_\phi(t)$ [rad]', 'Interpreter', 'latex');
title('Angle Error: $e_\phi(t) = \phi(t) - \hat{\phi}(t)$ (with $d\neq 0$)', 'Interpreter', 'latex');
yline(0, 'k--', 'LineWidth', 0.8);
grid on;

sgtitle('Exercise 1b: Gradient Method — Angle Estimation ($d(t) = 0.02\sin(2t)$)', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Plot 2: Parameter Estimates — d=0 vs d≠0 Comparison
figure('Name', 'Exercise 1b — Parameter Estimates Comparison', ...
       'NumberTitle', 'off', 'Position', [80, 100, 900, 600]);

subplot(2, 1, 1);
plot(t,   k_hat,   'b',  'LineWidth', 1.5); hold on;
plot(t_d, k_hat_d, 'r--','LineWidth', 1.5);
yline(k, 'k:', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\hat{k}(t)$', 'Interpreter', 'latex');
title('Parameter Estimate: $\hat{k}(t)$ — Comparison', 'Interpreter', 'latex');
legend('$\hat{k}(t)$, $d=0$', '$\hat{k}(t)$, $d\neq 0$', 'True $k$', ...
       'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t,   b_hat,   'b',  'LineWidth', 1.5); hold on;
plot(t_d, b_hat_d, 'r--','LineWidth', 1.5);
yline(b, 'k:', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$\hat{b}(t)$', 'Interpreter', 'latex');
title('Parameter Estimate: $\hat{b}(t)$ — Comparison', 'Interpreter', 'latex');
legend('$\hat{b}(t)$, $d=0$', '$\hat{b}(t)$, $d\neq 0$', 'True $b$', ...
       'Interpreter', 'latex', 'Location', 'best');
grid on;

sgtitle('Exercise 1b: Gradient Method — Parameter Estimates: $d=0$ vs $d(t)=0.02\sin(2t)$', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Plot 3: Parameter Errors — d=0 vs d≠0 Comparison
figure('Name', 'Exercise 1b — Parameter Errors Comparison', ...
       'NumberTitle', 'off', 'Position', [80, 100, 900, 600]);

subplot(2, 1, 1);
plot(t,   e_k,   'b',  'LineWidth', 1.5); hold on;
plot(t_d, e_k_d, 'r--','LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_k(t) = \hat{k}(t) - k$', 'Interpreter', 'latex');
title('Parameter Error $e_k(t)$: $d=0$ vs $d(t)\neq 0$', 'Interpreter', 'latex');
legend('$d = 0$', '$d(t) = 0.02\sin(2t)$', 'Interpreter', 'latex', 'Location', 'best');
yline(0, 'k--', 'LineWidth', 0.8);
grid on;

subplot(2, 1, 2);
plot(t,   e_b,   'b',  'LineWidth', 1.5); hold on;
plot(t_d, e_b_d, 'r--','LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$e_b(t) = \hat{b}(t) - b$', 'Interpreter', 'latex');
title('Parameter Error $e_b(t)$: $d=0$ vs $d(t)\neq 0$', 'Interpreter', 'latex');
legend('$d = 0$', '$d(t) = 0.02\sin(2t)$', 'Interpreter', 'latex', 'Location', 'best');
yline(0, 'k--', 'LineWidth', 0.8);
grid on;

sgtitle('Exercise 1b: Gradient Method — Parameter Errors: $d=0$ vs $d(t)=0.02\sin(2t)$', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

%% Summary
fprintf('===== Summary =====\n');
fprintf('Case d=0: k_hat=%.5f (|ek|=%.5f), b_hat=%.5f (|eb|=%.5f)\n', ...
        k_hat(end), abs(k - k_hat(end)), b_hat(end), abs(b - b_hat(end)));
fprintf('Case d≠0: k_hat=%.5f (|ek|=%.5f), b_hat=%.5f (|eb|=%.5f)\n', ...
        k_hat_d(end), abs(k - k_hat_d(end)), b_hat_d(end), abs(b - b_hat_d(end)));