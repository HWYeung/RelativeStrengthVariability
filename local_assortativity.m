function [Local_assort] = local_assortativity(W)
% Compute local assortativity for a 3D adjacency matrix W (R x R x N)
% Output: Local_assort (R x N) — node-level assortativity per subject

R = size(W, 1);
N = size(W, 3);

M = sum(W, [1 2]); % total edge weight per subject, 1 x 1 x N

% Node degrees per subject
j_i = sum(W, 1);   % 1 x R x N (row sum per column)
k_i = sum(W, 2);   % R x 1 x N (col sum per row)

% Broadcast degrees to match W shape (R x R x N)
J = repmat(j_i, [R, 1, 1]);
K = repmat(k_i, [1, R, 1]);

% Mean degree-related term
mean_deg = sum(0.5 * (J + K) .* W, [1 2]) ./ M;
mean_deg_sq = mean_deg.^2;

% Denominator
deg_sq = sum(0.5 * (J.^2 + K.^2) .* W, [1 2]) ./ M;
fixed_bottom = deg_sq - mean_deg_sq;

% Avoid divide-by-zero
fixed_bottom(fixed_bottom == 0) = eps;

% Numerator
numerator = sum(J .* K .* W, [1 2]) - j_i .* mean_deg_sq;

% Final local assortativity
Local_assort = squeeze(numerator ./ fixed_bottom);
end
