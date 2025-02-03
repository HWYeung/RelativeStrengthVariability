function [X_resid] = SequentialResid(X)
X_resid = X;
for xpos = 2:size(X,2)
    Mdl=fitglm(X_resid(:,1:xpos-1),X_resid(:,xpos));
    X_resid(:,xpos) = Mdl.Residuals.Raw;
end
X_resid = normalize(X_resid);