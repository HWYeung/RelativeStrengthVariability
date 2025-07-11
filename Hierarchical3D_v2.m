function [Rm,R] = Hierarchical3D_v2(W)
%Parameters:    W = Structural Connectome R x R x N
%
%Output:        Rm = averaged hierachical complexity
%               R = hierachical complexity for the corresponding
%                   node degree

n = size(W,1);

% this ensure that we are dealing with binary graph
W = W > 0; 


W2 = (sum(W,2)-zeros(size(W)));
NodeStrength = sum(W);
NodeStrengthcell = squeeze(num2cell(NodeStrength,[1 2]));
GroupCell = cellfun(@findgroups,NodeStrengthcell,'UniformOutput',false);
clear NodeStrength
W2 = sort(W2.* W);
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

clear Average
