function [t, x] = simulate_unknown(theta1, theta2, theta3, t_end, dt, u_func, x0)
%   Simulates the scalar nonlinear system whose structure 
%   is "unknown" to the estimator in Part 2.
%
%   System equation:  x_dot(t) = f(x(t), u(t))
%   where f(x,u) = -0.8 * x^3 * exp(x) + theta1*x 
%                + theta2*sin(x) + theta3*exp(-x) + u
%
%   Inputs:
%     theta1 - parameter on the x term       (in [-1.5, -0.5])
%     theta2 - parameter on the sin(x) term  (in [ 2.0,  4.0])
%     theta3 - parameter on the exp(-x) term (in [ 0.2,  0.8])
%     t_end  - simulation end time [s]
%     dt     - integration step [s]
%     u_func - input signal function handle @(t)
%     x0     - initial condition (scalar)
%
%   Outputs:
%     t - time vector [s]
%     x - state trajectory x(t)

% True dynamics
f = @(t, x) -0.8 * x^3 * exp(x) ...
            + theta1 * x ...
            + theta2 * sin(x) ...
            + theta3 * exp(-x) ...
            + u_func(t);

% Solver settings
t_span = 0 : dt : t_end;
options = odeset('MaxStep', dt, 'RelTol', 1e-6, 'AbsTol', 1e-9);

% Run ode45
[t, x] = ode45(f, t_span, x0, options);

end
