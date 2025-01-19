function [Rm,R] = Hierarchical3D(W,tiered,weighted, normalisation)
%Parameters:    W = Structural Connectome R x R x N
%               weighted = 1 for weighted average, 0 for unweighted
%               (default 0)
%               direction = 1 for arithmetic mean, 2 for harmonic mean
%Output:        NodeRelStrength = overall node strength relative to its
%                         neighbourhood

%If for arithmetic mean, the output will be inverted afterwards so 
% descend is actually correct here
n = size(W,1);
if nargin < 2
    tiered = 0;
end

if nargin < 3
    weighted = 0;
end


if nargin < 4
    normalisation = 0;
end



if normalisation == 1
    W = W./max(max(abs(W)));
end

W2 = (sum(W,2)-zeros(size(W)));
NodeStrength = squeeze(sum(W>0));
if tiered == 0
    NodeStrengthcell = num2cell(NodeStrength,[1 2]);
    GroupCell = cellfun(@findgroups,NodeStrengthcell,'UniformOutput',false);
else
    Quantiles = quantile(NodeStrength,3);
    DegTiers = 1;
    for tiers = 1:3
        DegTiers = DegTiers + (NodeStrength < Quantiles(tiers,:));
    end
    DegTierscell = num2cell(DegTiers,1);
    GroupCell = cellfun(@findgroups,DegTierscell,'UniformOutput',false);

end
clear NodeStrength
%W2 = W2;%./permute(W2,[2 1 3]);
if weighted == 0
    MeanCal = (W>0);
else
    MeanCal = W;
end
Average = squeeze(sum(MeanCal.*W2)./sum(MeanCal));
Averagecell = num2cell(Average,[1 2]);
countfun = @(y) splitapply(@sum,ones(1,n),y);
F1 = @(x,y) splitapply(@nanvar,x,y)/n;
F2 = @(x,y) nanmean(x(y>1));
takemean = cellfun(countfun,GroupCell,'UniformOutput',false);
R = cellfun(F1,Averagecell,GroupCell,'UniformOutput',false);
Rm = cellfun(F2,R,takemean)';

%StrengthVariability = std(NodeRelStrength)';
clear Average