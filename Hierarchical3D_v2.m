function [Rm,R] = Hierarchical3D_v2(W,weighted, normalisation)
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
    weighted = 0;
end


if nargin < 3
    normalisation = 0;
end

if normalisation == 1
    W = W./max(max(abs(W)));
end

W2 = (sum(W,2)-zeros(size(W)));
NodeStrength = sum(W);
NodeStrengthcell = squeeze(num2cell(NodeStrength,[1 2]));
GroupCell = cellfun(@findgroups,NodeStrengthcell,'UniformOutput',false);
clear NodeStrength
%W2 = W2;%./permute(W2,[2 1 3]);
if weighted == 0
    MeanCal = (W>0);
else
    MeanCal = W;
end
W2 = sort(W2.* MeanCal);
W2(W2 == 0) = nan;
W2cell = squeeze(num2cell(W2,[1 2]));
var2_fun = @(x) nanvar(x,[],2);
countfun = @(y) splitapply(@sum,ones(1,n),y);
F1 = @(x,y) splitapply(var2_fun,x,y)/n;
F2 = @(x) nanmean(x);
F3 = @(x,y) nanmean(x(y>1));
temptR = cellfun(F1,W2cell,GroupCell,'UniformOutput',false);
takemean = cellfun(countfun,GroupCell,'UniformOutput',false);
R = cellfun(F2,temptR,'UniformOutput',false);
Rm = cellfun(F3,R,takemean);

%StrengthVariability = std(NodeRelStrength)';
clear Average