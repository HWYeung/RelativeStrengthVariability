function [NodeRelStrength, StrengthVariability, windowed_StrVar] = NodeRelStrengthVariability(W,weighted, direction, normalisation,window_size)
%Parameters:    W = Structural Connectome R x R x N
%               weighted = 1 for weighted average, 0 for unweighted
%               (default 0)
%               direction = 1 for arithmetic mean, 2 for harmonic mean
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
%{
Quantiles = quantile(nodedegree,tiernum);
DegTiers = 1;
for tiers = 1:tiernum
    DegTiers = DegTiers + (nodedegree < Quantiles(tiers,:));
end
%}

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

%NodeRelStrengthcell = num2cell(NodeRelStrength,1);
%DegTierscell = num2cell(DegTiers,1);
%GroupCell = cellfun(@findgroups,DegTierscell,'UniformOutput',false);
%F1 = @(x,y) splitapply(@nanvar,x,y);
%F2 = @(x) nansum(x)./nansum(x>0);
%TieredVariability = cellfun(F1,NodeRelStrengthcell,GroupCell,'UniformOutput',false);
%TieredVariabilityM = cellfun(F2,TieredVariability');


clear Average
