function [NodeRelStrength, StrengthVariability, windowed_StrVar, hierarchical_variability] = ...
    NodeRelStrengthVariability(W, weighted, mean_type, normalisation, window_frac)

% Default arguments
if nargin < 2, weighted = 0; end
if nargin < 3, mean_type = 1; end
if nargin < 4, normalisation = 0; end
if nargin < 5, window_frac = 1; end

[R, ~, N] = size(W);
if normalisation == 1
    W = W ./ max(max(abs(W)));  % Normalize across all subjects
end

% Compute node relative strength
TotalStrength = sum(W, 2);           % Node strength to all others
W2 = TotalStrength - W;              % Neighborhood strength excluding self

% Normalize strength by neighborhood
denom = permute(W2, [2 1 3]);
W2 = W2 ./ (denom + (denom == 0));
W2(denom == 0) = 0;

valid_mask = (W > 0);
Average = sum(valid_mask .* W2, mean_type) ./ sum(valid_mask, mean_type);

% Compute relative strength
if mean_type == 1
    NodeRelStrength = 1 ./ squeeze(Average);
else
    NodeRelStrength = squeeze(Average);
end

% Compute variability
if weighted == 0
    StrengthVariability = std(NodeRelStrength)';
else
    total_strength = squeeze(sum(W, 1));
    weight = total_strength ./ sum(total_strength);
    mu = sum(NodeRelStrength .* weight, 2);
    StrengthVariability = sqrt(sum(weight .* (NodeRelStrength - mu).^2, 2));
end

% Hierarchical variability
window_size = floor(R * window_frac);
windowed_StrVar = [];
hierarchical_variability = [];

if window_size < R
    mean_strength = mean(squeeze(sum(W,1)), 2);
    [~, idx] = sort(mean_strength);
    sorted_Str = NodeRelStrength(idx,:);
    windowed_StrVar = zeros(N, R - window_size + 1);
    
    for i = 1:(R - window_size + 1)
        windowed_StrVar(:,i) = std(sorted_Str(i:i+window_size-1,:), 0, 1);
    end
    
    hierarchical_variability = mean(windowed_StrVar, 2);
end

