%%  Simulation and Modelling of Dynamic Systems
%  Lab 01 - Exercise 1

clear; clc; close all;

%% System Parameters
% Values chosen within the required ranges:
% m in [0.75, 1], k in [10, 15], c in [0.15, 0.35]
m = 0.9;    % mass [kg]
k = 12.0;   % spring constant [N/m]
c = 0.25;   % damping coefficient [N·s/m]

fprintf('===== System Parameters =====\n');
fprintf('m = %.2f kg\n', m);
fprintf('k = %.2f N/m\n', k);
fprintf('c = %.2f N·s/m\n\n', c);

%% State-Space Representation
% State variables: x1 = x(t) [position], x2 = x_dot(t) [velocity]
%
%  [x1_dot]   [  0    1  ] [x1]   [  0  ]
%  [x2_dot] = [-k/m  -c/m] [x2] + [ 1/m ] * u(t)
%
A = [0,    1;
    -k/m, -c/m];

B = [0;
     1/m];

C = [1, 0];
D = 0;

fprintf('===== State Matrix A =====\n'); disp(A)
fprintf('===== Input Matrix B =====\n'); disp(B)

%% Transfer Function
% Applying Laplace transform (zero initial conditions):
% G(s) = X(s)/U(s) = 1 / (m*s^2 + c*s + k)
num = [1];
den = [m, c, k];
G = tf(num, den);

fprintf('===== Transfer Function G(s) =====\n');
display(G)

%% Stability Check via Eigenvalues
poles = eig(A);
fprintf('System poles:\n'); disp(poles)
if all(real(poles) < 0)
    fprintf('System is stable (all poles have negative real part).\n\n');
else
    fprintf('WARNING: System may be unstable.\n\n');
end

%% Simulation Settings
dt = 1e-4;   % integration step (< 1e-3)
t_end = 20;  % simulation duration

%% Run Simulation
fprintf('Running ode45 simulation...\n');
[t, x1, x2] = simulate_system(m, k, c, t_end, dt);
fprintf('Total points: %d\n\n', length(t));

% Reconstruct input signal
u_out = 5 * sin(2.5 * t);

%% Plot: System States and Input
figure('Name', 'Exercise 1 — System Response', ...
       'NumberTitle', 'off', 'Position', [100, 100, 900, 650]);

% Position x(t)
subplot(3, 1, 1);
plot(t, x1, 'b', 'LineWidth', 1.5);
xlabel('Time t [s]', 'Interpreter', 'latex');
ylabel('x(t) [m]', 'Interpreter', 'latex');
title('State 1: Position x(t)', 'Interpreter', 'latex');;
grid on;

% Velocity x_dot(t)
subplot(3, 1, 2);
plot(t, x2, 'r', 'LineWidth', 1.5);
xlabel('Time t [s]', 'Interpreter', 'latex');
ylabel('$\dot{x}(t)$ [m/s]', 'Interpreter', 'latex');
title('State 2: Velocity $\dot{x}(t)$', 'Interpreter', 'latex');
grid on;

% Input u(t)
subplot(3, 1, 3);
plot(t, u_out,'g', 'LineWidth', 1.5);
xlabel('Time t [s]', 'Interpreter', 'latex');
ylabel('u(t) [N]', 'Interpreter', 'latex');
title('Input $u(t) = 5sin(2.5t)$', 'Interpreter', 'latex');
grid on;

sgtitle('Exercise 1: System Response', ...
        'FontSize', 13, 'FontWeight', 'bold');

%% Plot: Phase Portrait
figure('Name', 'Phase Portrait', 'NumberTitle', 'off', ...
       'Position', [100, 100, 600, 500]);

plot(x1, x2, 'm', 'LineWidth', 1.2);
xlabel('$x(t)$ [m]', 'Interpreter', 'latex');
ylabel('$\dot{x}(t)$ [m/s]', 'Interpreter', 'latex');
title('Phase Portrait: Velocity vs Position');
grid on;

%% Summary Statistics
fprintf('===== Simulation Summary =====\n');
fprintf('Peak position: %.4f m\n',   max(abs(x1)));
fprintf('Peak velocity: %.4f m/s\n', max(abs(x2)));

% Save workspace data
save('simulation_data.mat', 't', 'x1', 'x2', 'u_out', 'dt', 'm', 'k', 'c');
fprintf('\nSimulation data saved.\n\n');