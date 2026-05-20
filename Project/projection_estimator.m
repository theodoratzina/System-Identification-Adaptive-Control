function [A_hat, b_hat, x1_hat, x2_hat] = projection_estimator(t, x1, x2, u, K, ...
    gamma_A, gamma_b, A0, b0, A_lower, A_upper, b_lower, b_upper)
%   Projection-based Lyapunov estimator for the 2-D nonlinear MIMO system
%   under unknown bounded bias disturbance.
%
%   System equation:  x_dot(t) = A*x(t) + G(x(t))*u(t) + w(t)
%   where G(x) = [b1 + sin(x1) + x2^2;
%                 b2 / (1 + x1^2)    ]
%   and  ||w(t)|| <= w_bar (w_bar > 0, unknown)
%
%   Inputs:
%     t       - time vector            [N x 1]
%     x1      - measured state x1(t)   [N x 1]
%     x2      - measured state x2(t)   [N x 1]
%     u       - control input          [N x 1]
%     K       - stabilizing gains      [2 x 1]
%     gamma_A - adaptation gains for A [4 x 1]
%     gamma_b - adaptation gains for b [2 x 1]
%     A0      - initial estimate of A  [2 x 2]
%     b0      - initial estimate of b  [2 x 1]
%     A_lower - lower bounds for A     [4 x 1]
%     A_upper - upper bounds for A     [4 x 1]
%     b_lower - lower bounds for b     [2 x 1]
%     b_upper - upper bounds for b     [2 x 1]
%
%   Outputs:
%     A_hat  - parameter history of A     [N x 4]
%     b_hat  - parameter history of b     [N x 2]
%     x1_hat - estimated state trajectory [N x 1]
%     x2_hat - estimated state trajectory [N x 1]

N = length(t);

A_hat = zeros(N, 4);
b_hat = zeros(N, 2);
x1_hat = zeros(N, 1);
x2_hat = zeros(N, 1);

% Initialise (must start strictly inside the box w)
A_hat(1, :) = [A0(1,1), A0(1,2), A0(2,1), A0(2,2)];
b_hat(1, :) = b0(:)';
x1_hat(1) = x1(1);
x2_hat(1) = x2(1);

for i = 1 : N - 1
    dt_i = t(i+1) - t(i);

    % Current parameter estimates as column vectors
    A_cur = A_hat(i, :)';
    b_cur = b_hat(i, :)';

    % Current true states and input
    x1_i = x1(i);
    x2_i = x2(i);
    u_i = u(i);

    % State tracking errors (Lyapunov driving signals)
    e1 = x1_i - x1_hat(i);
    e2 = x2_i - x2_hat(i);

    % Series-parallel model dynamics
    x1_hat_dot = A_cur(1)*x1_i + A_cur(2)*x2_i + b_cur(1)*u_i ...
                 + (sin(x1_i) + x2_i^2)*u_i + K(1)*e1;
    x2_hat_dot = A_cur(3)*x1_i + A_cur(4)*x2_i ...
                 + b_cur(2)*u_i / (1 + x1_i^2) + K(2)*e2;

    % Euler integration of parallel-model states
    x1_hat(i+1) = x1_hat(i) + dt_i * x1_hat_dot;
    x2_hat(i+1) = x2_hat(i) + dt_i * x2_hat_dot;

    % Standard (unprojected) Lyapunov update directions
    tau_A = [gamma_A(1) * e1 * x1_i;
             gamma_A(2) * e1 * x2_i;
             gamma_A(3) * e2 * x1_i;
             gamma_A(4) * e2 * x2_i];

    tau_b = [gamma_b(1) * e1 * u_i;
             gamma_b(2) * e2 * u_i / (1 + x1_i^2)];

    % Apply box-projection element-wise
    tau_A_proj = projection_operator(A_cur, tau_A, A_lower, A_upper);
    tau_b_proj = projection_operator(b_cur, tau_b, b_lower, b_upper);

    % Euler step on projected update directions
    A_hat(i+1, :) = A_cur' + dt_i * tau_A_proj';
    b_hat(i+1, :) = b_cur' + dt_i * tau_b_proj';

    % Numerical safety net: clamp in case discrete Euler step crosses bound
    A_hat(i+1, :) = max(min(A_hat(i+1, :), A_upper'), A_lower');
    b_hat(i+1, :) = max(min(b_hat(i+1, :), b_upper'), b_lower');
end

end


function tau_proj = projection_operator(theta_hat, tau, theta_min, theta_max)
%   Box-projection operator for adaptive parameter estimation.
%
%   The operator satisfies the key property:
%       (theta_hat - theta*)' * Gamma^(-1) * (Proj(tau) - tau) <= 0
%
%   Inputs:
%     theta_hat - current parameter estimate [n x 1]
%     tau       - proposed update direction  [n x 1]
%     theta_min - lower bound, element-wise  [n x 1]
%     theta_max - upper bound, element-wise  [n x 1]
%
%   Output:
%     tau_proj - projected update direction [n x 1]

tau_proj = tau;

% Upper bound and update direction points further out → kill it
mask_upper = (theta_hat >= theta_max) & (tau > 0);
tau_proj(mask_upper) = 0;

% Lower bound and update direction points further out → kill it
mask_lower = (theta_hat <= theta_min) & (tau < 0);
tau_proj(mask_lower) = 0;

end