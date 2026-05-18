function [t, x1, x2] = simulate_mimo(A, b1, b2, t_end, dt, u_func, w_func)
%   Simulates the 2-D nonlinear MIMO system.
%
%   System equation:  x_dot(t) = A*x(t) + G(x(t))*u(t) + w(t)
%   where G(x) = [b1 + sin(x1) + x2^2;
%                 b2 / (1 + x1^2)    ]
%
%   Inputs:
%     A      - 2x2 system matrix
%     b1     - input nonlinearity parameter (channel 1)
%     b2     - input nonlinearity parameter (channel 2)
%     t_end  - simulation end time [s]
%     dt     - integration step [s]
%     u_func - input signal function handle @(t)
%     w_func - bias / disturbance function handle @(t), returns [2x1] vector
%              (use @(t) [0; 0] for the disturbance-free case)
%
%   Outputs:
%     t  - time vector [s]
%     x1 - state 1: x1(t)
%     x2 - state 2: x2(t)
%
%   State variables:
%       x1 = x1(t)
%       x2 = x2(t)
%
%   State-space form:
%       x1_dot = a11*x1 + a12*x2 + (b1 + sin(x1) + x2^2)*u(t) + w1(t)
%       x2_dot = a21*x1 + a22*x2 + (b2/(1 + x1^2))*u(t) + w2(t)

% Nonlinear input gain G(x)
G = @(x) [b1 + sin(x(1)) + x(2)^2;
          b2 / (1 + x(1)^2)];

% ODE-Solver function
odefun = @(t, x) A*x + G(x)*u_func(t) + w_func(t);

% Initial conditions (zero)
x0 = [0; 0];

% Solver options
t_span = 0 : dt : t_end;
options = odeset('MaxStep', dt, 'RelTol', 1e-6, 'AbsTol', 1e-9);

% Run ode45
[t, X] = ode45(odefun, t_span, x0, options);

% Extract states
x1 = X(:, 1);
x2 = X(:, 2);

end
