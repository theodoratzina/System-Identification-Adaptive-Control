function [t, x1, x2] = simulate_system(m, k, c, t_end, dt)
%   Simulates the system using ode45.
%
%   Inputs:
%       m     - mass [kg]
%       k     - spring constant [N/m]
%       c     - damping coefficient [N·s/m]
%       t_end - simulation end time [s]
%       dt    - integration step size [s]
%
%   Outputs:
%       t     - time vector [s]
%       x1    - position x(t) [m]
%       x2    - velocity x_dot(t) [m/s]

% State variables: x1 = x(t), x2 = x_dot(t)
% x1_dot =  x2
% x2_dot = -k/m * x1 - c/m * x2 + 1/m * u(t)
A = [0,    1;
    -k/m, -c/m];

B = [0;
     1/m];

% Input signal
u = @(t) 5 * sin(2.5 * t);

% ODE-Solver function
odefun = @(t, x) A*x + B*u(t);

% Initial conditions (zero)
x0 = [0; 0];

% Solver options
t_span   = 0 : dt : t_end;
options = odeset('MaxStep', dt, 'RelTol', 1e-6, 'AbsTol', 1e-9);

% Run ode45
[t, X] = ode45(odefun, t_span, x0, options);

% Extract states
x1 = X(:, 1);   % position  [m]
x2 = X(:, 2);   % velocity  [m/s]

end