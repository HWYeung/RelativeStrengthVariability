function [NodeRelStrength, StrengthVariability, windowed_StrVar, hierarchical_variability] = NodeRelStrengthVariability(W,weighted, direction, normalisation,window_size)
%Parameters:    W = Structural Connectome R x R x N
%               weighted = 1 for weighted average, 0 for unweighted
%               (default 0)
%               direction = 1 for arithmetic mean, 2 for harmonic mean
%               normalisation = 1 for normalisation to [0, 1], 0 for no normalisation (default 0)
%               window_size = n<= number of nodes in graph, window size for calculating hierarchical variability measure
%Output:        NodeRelStrength = overall node strength relative to its
%                         neighbourhood

%If for arithmetic mean, the output will be inverted afterwards so 
% descend is actually correct here
nodedegree = squeeze(sum(W>0));

if nargin < 2
    weighted = 0;
end

if nargin < 3
    direction = 1;
end

if nargin < 4
    normalisation = 0;
end

if nargin < 5
    window_size = size(w,1);
end

if normalisation == 1
    W = W./max(max(abs(W)));
end

W2 = (sum(W,2)-W);
W2 = W2./permute(W2,[2 1 3]);
MeanCal = (W>0);

Average = sum(MeanCal.*W2,direction)./sum(MeanCal,direction);
clear W2 MeanCal
if direction == 1
    NodeRelStrength = 1./squeeze(Average);
elseif direction == 2
    NodeRelStrength = squeeze(Average);
end
clear Average
if weighted == 0
    StrengthVariability = std(NodeRelStrength)';
else
    Strength = squeeze(sum(W,1));
    Weightings = Strength./sum(Strength);
    WeightedMean = sum(NodeRelStrength.*Weightings);
    SquaredNodeRel = (NodeRelStrength - WeightedMean).^2;
    StrengthVariability = sum(Weightings.*SquaredNodeRel).^(0.5);
end
if window_size == size(W,1)
    return
end
MeanNodeStrength = mean(squeeze(sum(W,1)),2);
[~,NodeStrengthOrder] = sort(MeanNodeStrength);
clear MeanNodeStrength
NodeRelStrength_sorted = NodeRelStrength(NodeStrengthOrder,:);
Window_size = floor(size(W,1)*window_size);
windowed_StrVar = zeros(size(W,3), size(W,1) - Window_size + 1);
for i = 1:size(windowed_StrVar,2)
    windowed_StrVar(:,i) = std(NodeRelStrength_sorted(i:i+Window_size-1,:));
end

hierarchical_variability = mean(windowed_StrVar, 2)

clear Average
