%% Simulation and Modelling of Dynamic Systems
% Final Project — Part 2
% Model Structure Selection and Generalisation for an Unknown Nonlinear System

clear; clc; close all;

%% True System Parameters
theta1_true = -1.0;
theta2_true = 3.0;
theta3_true = 0.5;

fprintf('===== True System (assumed unknown to the estimator) =====\n');
fprintf('x_dot = -0.8*x^3*exp(x) + (%.2f)*x + (%.2f)*sin(x) + (%.2f)*exp(-x) + u\n\n', ...
        theta1_true, theta2_true, theta3_true);

%% Simulation Settings
dt = 1e-3;
t_end = 150;   % training/validation/test horizon [s]

%% Dataset Generation
% Training input: rich tri-tone for strong PE
u_train_func = @(t) 2.0*sin(0.7*t) + 1.5*sin(1.5*t) + 1.0*cos(0.3*t) + ...
                    1.0*sin(2.3*t) + 0.7*cos(3.1*t);
x0_train = 0;

% Validation input: different frequencies/amplitudes
u_val_func = @(t) 1.5*sin(0.5*t) + 1.0*cos(1.2*t) + 1.5*sin(2.0*t) + ...
                  0.8*cos(2.8*t) + 0.6*sin(3.7*t);
x0_val = 0.3;

% Test input: completely different + different IC for generalisation
u_test_func = @(t) 2.5*cos(0.4*t) + 0.8*sin(1.8*t) + 1.0*cos(0.9*t) + ...
                   0.8*sin(2.5*t) + 0.5*cos(3.4*t);
x0_test = -0.5;

fprintf('===== Generating Datasets =====\n');
[t_train, x_train] = simulate_unknown(theta1_true, theta2_true, theta3_true, ...
                                      t_end, dt, u_train_func, x0_train);
u_train = u_train_func(t_train);

[t_val, x_val] = simulate_unknown(theta1_true, theta2_true, theta3_true, ...
                                  t_end, dt, u_val_func, x0_val);
u_val = u_val_func(t_val);

[t_test, x_test] = simulate_unknown(theta1_true, theta2_true, theta3_true, ...
                                    t_end, dt, u_test_func, x0_test);
u_test = u_test_func(t_test);

fprintf('State ranges:\n');
fprintf('  Training:   x in [%+.2f, %+.2f]\n', min(x_train), max(x_train));
fprintf('  Validation: x in [%+.2f, %+.2f]\n', min(x_val), max(x_val));
fprintf('  Test:       x in [%+.2f, %+.2f]\n\n', min(x_test), max(x_test));

%% Figure 1: Dataset Overview
figure('Name', 'Exercise 2 - Dataset Overview', ...
       'NumberTitle', 'off', 'Position', [80, 100, 1200, 700]);

subplot(2, 3, 1);
plot(t_train, x_train, 'b', 'LineWidth', 1.2);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$x_{train}(t)$', 'Interpreter', 'latex');
title('Training Data', 'Interpreter', 'latex');
grid on;

subplot(2, 3, 2);
plot(t_val, x_val, 'g', 'LineWidth', 1.2);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$x_{val}(t)$', 'Interpreter', 'latex');
title('Validation Data', 'Interpreter', 'latex');
grid on;

subplot(2, 3, 3);
plot(t_test, x_test, 'r', 'LineWidth', 1.2);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$x_{test}(t)$', 'Interpreter', 'latex');
title('Test Data', 'Interpreter', 'latex');
grid on;

subplot(2, 3, 4);
plot(t_train, u_train, 'b', 'LineWidth', 1.2);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$u_{train}(t)$', 'Interpreter', 'latex');
title('Training Input', 'Interpreter', 'latex');
grid on;

subplot(2, 3, 5);
plot(t_val, u_val, 'g', 'LineWidth', 1.2);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$u_{val}(t)$', 'Interpreter', 'latex');
title('Validation Input', 'Interpreter', 'latex');
grid on;

subplot(2, 3, 6);
plot(t_test, u_test, 'r', 'LineWidth', 1.2);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$u_{test}(t)$', 'Interpreter', 'latex');
title('Test Input', 'Interpreter', 'latex');
grid on;

sgtitle('Exercise 2: Training / Validation / Test Datasets', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

%% Exercise 2a — Structure Selection & Cross-Validation
% Stabilising gain (same for all candidates)
K = 10;

% Each model: name, basis function handle phi(x,u), initial parameter
% estimate, and per-parameter adaptation gain vector. Gains are individually
% scaled to compensate for the very different amplitude ranges of the bases
M1 = struct( ...
    'name',   'M1: Linear', ...
    'phi',    @(x,u) [x; u], ...
    'theta0', zeros(2, 1), ...
    'gamma',  [150; 150] );

M2 = struct( ...
    'name',   'M2: Polynomial', ...
    'phi',    @(x,u) [x; x^2; x^3; u], ...
    'theta0', zeros(4, 1), ...
    'gamma',  [150; 100; 60; 150] );

M3 = struct( ...
    'name',   'M3: Correct structure', ...
    'phi',    @(x,u) [x; sin(x); exp(-x); x^3*exp(x); u], ...
    'theta0', zeros(5, 1), ...
    'gamma',  [150; 150; 30; 30; 150] );

M4 = struct( ...
    'name',   'M4: Overparameterised', ...
    'phi',    @(x,u) [x; x^2; x^3; sin(x); cos(x); ...
                      exp(-x); exp(x); x^3*exp(x); u], ...
    'theta0', zeros(9, 1), ...
    'gamma',  [150; 100; 60; 150; 150; 30; 100; 30; 150] );

models = {M1, M2, M3, M4};
N_models = length(models);

fprintf('===== Candidate Models =====\n');
for j = 1:N_models
    fprintf('  %s (M = %d parameters)\n', models{j}.name, length(models{j}.theta0));
end
fprintf('\n');

%% Train Each Model + Cross-Validate + Test
results = struct( ...
    'theta_final', cell(N_models, 1), ...
    'x_hat_train', cell(N_models, 1), ...
    'x_sim_val',   cell(N_models, 1), ...
    'x_sim_test',  cell(N_models, 1), ...
    'K_train',     cell(N_models, 1), ...
    'K_val',       cell(N_models, 1), ...
    'K_test',      cell(N_models, 1) );

fprintf('===== Training and Cross-Validation =====\n');
fprintf('%-28s %-4s %-12s %-12s %-12s\n', ...
        'Model', 'M', 'K_train', 'K_val', 'K_test');
fprintf('%s\n', repmat('-', 1, 72));

for j = 1:N_models
    Mj = models{j};

    % 1. Train (adaptive series-parallel estimator)
    [theta_h, x_hat_tr] = basis_estimator(t_train, x_train, u_train, ...
                                          Mj.phi, Mj.gamma, K, Mj.theta0);
    results(j).theta_final = theta_h(end, :)';
    results(j).x_hat_train = x_hat_tr;
    results(j).K_train = trapz(t_train, (x_train - x_hat_tr).^2);

    % 2. Validate (pure-parallel simulation, frozen params)
    x_sim_val = pure_simulate(t_val, u_val, Mj.phi, results(j).theta_final, x_val(1));
    results(j).x_sim_val = x_sim_val;
    results(j).K_val = trapz(t_val, (x_val - x_sim_val).^2);

    % 3. Test (generalisation)
    x_sim_test = pure_simulate(t_test, u_test, Mj.phi, results(j).theta_final, x_test(1));
    results(j).x_sim_test = x_sim_test;
    results(j).K_test = trapz(t_test, (x_test - x_sim_test).^2);

    fprintf('%-28s %-4d %-12.4f %-12.4f %-12.4f\n', ...
            Mj.name, length(Mj.theta0), ...
            results(j).K_train, results(j).K_val, results(j).K_test);
end
fprintf('\n');

%% Exercise 2b - Final Model Selection & Generalisation Test
% Minimisation of modelling error K vs complexity
K_train_arr = [results.K_train];
K_val_arr = [results.K_val];
K_test_arr = [results.K_test];
M_arr = arrayfun(@(j) length(models{j}.theta0), 1:N_models);

% Pick lowest validation error first
[K_val_min, best_idx] = min(K_val_arr);
fprintf('===== Model Selection =====\n');
fprintf('Lowest validation error: %s with K_val = %.4f\n', ...
         models{best_idx}.name, K_val_min);

% Parsimony principle: prefer simpler model if K_val is within 10% of the minimum
tolerance = 1.10;
mask = (K_val_arr <= K_val_min * tolerance) & (M_arr < M_arr(best_idx));
candidate_simpler = find(mask);

if ~isempty(candidate_simpler)
    [~, k] = min(M_arr(candidate_simpler));
    final_idx = candidate_simpler(k);
    fprintf('Parsimony: model %s achieves K_val = %.4f with fewer parameters (M = %d).\n', ...
            models{final_idx}.name, K_val_arr(final_idx), M_arr(final_idx));
    fprintf('Selecting it as the FINAL model.\n');
else
    final_idx = best_idx;
    fprintf('No simpler model is within %.0f%% of the minimum K_val. Keeping %s.\n', ...
            (tolerance-1)*100, models{final_idx}.name);
end
fprintf('\n');

%% Figure 2: Cross-Validation Performance
figure('Name', 'Exercise 2 - Cross-Validation Performance', ...
       'NumberTitle', 'off', 'Position', [80, 80, 900, 600]);

bar_data = [K_train_arr; K_val_arr; K_test_arr]';
b_h = bar(bar_data);
set(gca, 'YScale', 'log');
set(gca, 'XTickLabel', arrayfun(@(j) models{j}.name, 1:N_models, ...
        'UniformOutput', false));
xtickangle(15);
ylabel('Modelling error $K = \int(x-\hat{x})^2\,dt$', 'Interpreter', 'latex');
title('Exercise 2: Cross-Validation Performance (log scale)', 'Interpreter', 'latex');
legend('$K_{train}$', '$K_{val}$', '$K_{test}$', ...
       'Interpreter', 'latex', 'Location', 'best');
grid on;

% Annotate each model with its parameter count
y_top = max(K_test_arr) * 3;
for j = 1:N_models
    text(j, y_top, sprintf('M = %d', M_arr(j)), ...
         'HorizontalAlignment', 'center', 'FontSize', 10, ...
         'FontWeight', 'bold');
end
ylim([min(K_train_arr)/2, y_top*2]);

%% Figure 3: Validation Tracking per Candidate Model
figure('Name', 'Exercise 2 - Validation Tracking', ...
       'NumberTitle', 'off', 'Position', [100, 60, 1100, 700]);

for j = 1:N_models
    subplot(2, 2, j);
    plot(t_val, x_val, 'b', 'LineWidth', 1.5); hold on;
    plot(t_val, results(j).x_sim_val, 'r--', 'LineWidth', 1.5);
    xlabel('Time $t$ [s]', 'Interpreter', 'latex');
    ylabel('$x(t)$', 'Interpreter', 'latex');
    title(sprintf('%s ($K_{val} = %.3f$)', models{j}.name, results(j).K_val), ...
          'Interpreter', 'latex');
    legend('True $x$', 'Model $\hat{x}$', ...
           'Interpreter', 'latex', 'Location', 'best');
    grid on;
end

sgtitle('Exercise 2: Pure-Simulation Tracking per Candidate Model', ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

%% Figure 4: Final Model Estimated Parameters
final_theta = results(final_idx).theta_final;
final_M = length(final_theta);

figure('Name', 'Exercise 2 - Final Model Parameters', ...
       'NumberTitle', 'off', 'Position', [120, 40, 800, 500]);

bar(1:final_M, final_theta, 'FaceColor', [0.2, 0.4, 0.8]);
xlabel('Basis function index $i$', 'Interpreter', 'latex');
ylabel('$\hat{\theta}_i$', 'Interpreter', 'latex');
title(sprintf('Exercise 2: Final Estimated Parameters - %s', models{final_idx}.name), ...
      'Interpreter', 'latex');
grid on;
xticks(1:final_M);

% Display values above each bar
for i = 1:final_M
    offset = sign(final_theta(i)) * 0.05 * max(abs(final_theta));
    if offset == 0; offset = 0.05 * max(abs(final_theta)); end
    text(i, final_theta(i) + offset, sprintf('%.3f', final_theta(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
end

%% Figure 5: Generalisation Test of Final Model
figure('Name', 'Exercise 2 - Generalisation Test', ...
       'NumberTitle', 'off', 'Position', [140, 20, 1000, 700]);

subplot(2, 1, 1);
plot(t_test, x_test, 'b', 'LineWidth', 1.5); hold on;
plot(t_test, results(final_idx).x_sim_test, 'r--', 'LineWidth', 1.5);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('$x(t)$', 'Interpreter', 'latex');
title(sprintf('Final model %s on TEST data ($K_{test} = %.3f$)', ...
              models{final_idx}.name, results(final_idx).K_test), ...
      'Interpreter', 'latex');
legend('True $x_{test}$', 'Model $\hat{x}_{test}$', ...
       'Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2, 1, 2);
plot(t_test, x_test - results(final_idx).x_sim_test, 'm', 'LineWidth', 1.5);
yline(0, 'k--', 'LineWidth', 0.8);
xlabel('Time $t$ [s]', 'Interpreter', 'latex');
ylabel('Generalisation error $e(t)$', 'Interpreter', 'latex');
title('Tracking Error on Test Data', 'Interpreter', 'latex');
grid on;

sgtitle(sprintf('Exercise 2: Generalisation Test: Final Model = %s', models{final_idx}.name), ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');

%% Summary
fprintf('===== Final Model Summary =====\n');
fprintf('Selected model: %s (M = %d parameters)\n', ...
        models{final_idx}.name, M_arr(final_idx));
fprintf('  K_train = %.4f\n', results(final_idx).K_train);
fprintf('  K_val   = %.4f\n', results(final_idx).K_val);
fprintf('  K_test  = %.4f\n', results(final_idx).K_test);

gen_ratio = results(final_idx).K_test / results(final_idx).K_val;
fprintf('  K_test / K_val = %.3f', gen_ratio);
if gen_ratio < 1.5
    fprintf('  (good generalisation)\n');
elseif gen_ratio < 3
    fprintf('  (acceptable generalisation)\n');
else
    fprintf('  (poor generalisation - possible overfitting)\n');
end
fprintf('\nFinal estimated parameters:\n');
disp(results(final_idx).theta_final');

fprintf('===== Part 2 Complete =====\n');

%% Local function: pure parallel simulation
function x_sim = pure_simulate(t, u_vec, phi, theta, x0)
%   Pure parallel simulation of candidate model with frozen parameters:
%       d/dt x_sim = theta^T * phi(x_sim, u)

t_end_local = t(end);
u_interp = @(tt) interp1(t, u_vec, tt, 'linear', 'extrap');

odefun = @(tt, xx) theta(:)' * phi(xx, u_interp(tt));

% Stop the solver if the trajectory explodes
options = odeset('Events', @explode_event, ...
                 'RelTol', 1e-6, 'AbsTol', 1e-9);

try
    sol = ode45(odefun, [0, t_end_local], x0, options);
    x_sim = deval(sol, t)';
catch
    % ode45 itself crashed - return saturated trajectory
    x_sim = sign(x0) * 100 * ones(size(t));
    if x0 == 0
        x_sim = 100 * ones(size(t));
    end
    return;
end

% If ode45 stopped early due to event, pad remaining samples with
% the saturation value
if length(x_sim) < length(t)
    last_val = x_sim(end);
    x_sim(end+1:length(t)) = sign(last_val) * 100;
end

% Clamp any NaN/Inf entries
bad = ~isfinite(x_sim);
if any(bad)
    x_sim(bad) = 100;
end
end


function [value, isterminal, direction] = explode_event(~, x)
value = 100 - abs(x);   % trigger when |x| reaches 100
isterminal = 1;         % stop integration
direction  = -1;        % only when decreasing through zero
end