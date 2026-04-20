function [t, phi, phi_dot, phi_ddot] = simulate_vehicle(J, k, b, t_end, dt, d_func)
%   Simulates the linearised single-axis orientation dynamics of an aerial vehicle.
%
%   System equation:  J * phi_ddot(t) = -k * phi_dot(t) + b * u(t) + d(t)
%
%   Inputs:
%     J      - moment of inertia [kg·m²]
%     k      - damping coefficient
%     b      - input gain constant
%     t_end  - simulation end time [s]
%     dt     - integration step [s]
%     d_func - disturbance function handle @(t)
%
%   Outputs:
%     t        - time vector [s]
%     phi      - orientation angle     [rad]
%     phi_dot  - angular velocity      [rad/s]
%     phi_ddot - angular acceleration  [rad/s²]
%
%   State variables:
%       x1 = phi(t)      [rad]    orientation angle
%       x2 = phi_dot(t)  [rad/s]  angular velocity
%
%   State-space form:
%       x1_dot =  x2
%       x2_dot = (-k*x2 + b*u(t) + d(t)) / J

% Input signal
u_func = @(t) 0.25 * sin(0.5 * pi * t);

% ODE-Solver function
odefun = @(t, x) [ x(2);
                   (-k * x(2) + b * u_func(t) + d_func(t)) / J ];

% Initial conditions (zero)
x0 = [0; 0];

% Solver options
t_span = 0 : dt : t_end;
options = odeset('MaxStep', dt, 'RelTol', 1e-6, 'AbsTol', 1e-9);

% Run ode45
[t, X] = ode45(odefun, t_span, x0, options);

% Extract states
phi = X(:, 1);      % orientation angle
phi_dot = X(:, 2);  % angular velocity

% Recompute acceleration analytically (avoiding numerical differentiation noise)
u = 0.25 * sin(0.5 * pi * t);
d = arrayfun(d_func, t);
phi_ddot = (-k * phi_dot + b * u + d) / J;

end
