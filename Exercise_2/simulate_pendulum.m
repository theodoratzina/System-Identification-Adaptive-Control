function [t, theta, theta_dot] = simulate_pendulum(g, l, c, t_end, dt)
%   Simulates the nonlinear simple pendulum.
%
%   System equation:  theta_ddot(t) = -(g/l)*sin(theta(t)) - c*theta_dot(t) + u(t)
%
%   Inputs:
%     g     - gravitational acceleration [m/s²]
%     l     - pendulum length [m]
%     c     - damping coefficient
%     t_end - simulation end time [s]
%     dt    - integration step [s]
%
%   Outputs:
%     t         - time vector [s]
%     theta     - pendulum angle [rad]
%     theta_dot - angular velocity [rad/s]
%
%   State variables:
%       x1 = theta(t)      [rad]    pendulum angle
%       x2 = theta_dot(t)  [rad/s]  angular velocity
%
%   State-space form:
%       x1_dot =  x2
%       x2_dot = -(g/l)*sin(x1) - c*x2 + u(t)

% Input signal
u_func = @(t) 0.5 * sin(t);

% ODE-Solver function
odefun = @(t, x) [ x(2);
                  -(g / l) * sin(x(1)) - c * x(2) + u_func(t) ];

% Small non-zero initial angle to excite the dynamics
x0 = [0.2; 0];

% Solver options
t_span = 0 : dt : t_end;
options = odeset('MaxStep', dt, 'RelTol', 1e-6, 'AbsTol', 1e-9);

% Run ode45
[t, X] = ode45(odefun, t_span, x0, options);

% Extract states
theta = X(:, 1);      % pendulum angle
theta_dot = X(:, 2);  % angular velocity

end
