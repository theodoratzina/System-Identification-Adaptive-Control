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
