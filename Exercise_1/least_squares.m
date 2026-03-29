function [m_hat, c_hat, k_hat] = least_squares(x1, x2, x2_dot, u)
%   Estimates system parameters via the Least Squares method.
%
%   Inputs:
%       x1    - sampled position x(t)          [N x 1]
%       x2    - sampled velocity x_dot(t)      [N x 1]
%       x2_dot- sampled acceleration x_ddot(t) [N x 1]
%       u     - sampled input u(t)             [N x 1]
%
%   Outputs:
%       m_hat - estimated mass [kg]
%       c_hat - estimated damping coefficient [N·s/m]
%       k_hat - estimated spring constant [N/m]

% Build regressor matrix Phi (N x 3)
Phi = [x2_dot, x2, x1];

% Output vector Y
Y = u;

% Least Squares solution
% theta = (Phi' * Phi)^-1 * Phi' * Y
theta_hat = (Phi' * Phi) \ (Phi' * Y);

% Extract individual parameter estimates
m_hat = theta_hat(1);
c_hat = theta_hat(2);
k_hat = theta_hat(3);

end