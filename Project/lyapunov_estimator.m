function [A_hat, b_hat, x1_hat, x2_hat] = lyapunov_estimator( ...
    t, x1, x2, u, K, gamma_A, gamma_b, A0, b0)
%   Implements the Lyapunov-based method for online estimation of the system
%   matrix A and the input nonlinearity parameters b1, b2.
%
%   System equation:  x_dot(t) = A*x(t) + G(x(t))*u(t)
%   where G(x) = [b1 + sin(x1) + x2^2;
%                 b2 / (1 + x1^2)    ]
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
%
%   Outputs:
%     A_hat  - parameter history of A     [N x 4]
%     b_hat  - parameter history of b     [N x 2]
%     x1_hat - estimated state trajectory [N x 1]
%     x2_hat - estimated state trajectory [N x 1]
%
%   The known nonlinear term (sin(x1)+x2^2)*u is fed forward directly — it has
%   no unknown parameter, so it is not subject to estimation.

N = length(t);

A_hat = zeros(N, 4);
b_hat = zeros(N, 2);
x1_hat = zeros(N, 1);
x2_hat = zeros(N, 1);

% Initialise — align parallel model with true initial states
A_hat(1, :) = [A0(1,1), A0(1,2), A0(2,1), A0(2,2)];
b_hat(1, :) = b0(:)';
x1_hat(1) = x1(1);
x2_hat(1) = x2(1);

for i = 1 : N - 1
    dt_i = t(i+1) - t(i);

    % Current parameter estimates
    a11 = A_hat(i, 1);
    a12 = A_hat(i, 2);
    a21 = A_hat(i, 3);
    a22 = A_hat(i, 4);
    b1 = b_hat(i, 1);
    b2 = b_hat(i, 2);

    % Current true states and input
    x1_i = x1(i);
    x2_i = x2(i);
    u_i = u(i);

    % State tracking errors (Lyapunov driving signals)
    e1 = x1_i - x1_hat(i);
    e2 = x2_i - x2_hat(i);

    % Series-parallel model dynamics
    % Regressor uses measured true states (not estimated ones)
    x1_hat_dot = a11*x1_i + a12*x2_i + b1*u_i + (sin(x1_i) + x2_i^2)*u_i + K(1)*e1;
    x2_hat_dot = a21*x1_i + a22*x2_i + b2*u_i / (1 + x1_i^2) + K(2)*e2;

    % Euler integration: advance parallel-model states
    x1_hat(i+1) = x1_hat(i) + dt_i * x1_hat_dot;
    x2_hat(i+1) = x2_hat(i) + dt_i * x2_hat_dot;

    % Lyapunov-based parameter update laws
    A_hat(i+1, 1) = a11 + dt_i * gamma_A(1) * e1 * x1_i;
    A_hat(i+1, 2) = a12 + dt_i * gamma_A(2) * e1 * x2_i;
    A_hat(i+1, 3) = a21 + dt_i * gamma_A(3) * e2 * x1_i;
    A_hat(i+1, 4) = a22 + dt_i * gamma_A(4) * e2 * x2_i;
    b_hat(i+1, 1) = b1  + dt_i * gamma_b(1) * e1 * u_i;
    b_hat(i+1, 2) = b2  + dt_i * gamma_b(2) * e2 * u_i / (1 + x1_i^2);
end

end
