function [lowApprox, PCA_score] = LRAMissing(X,rank)
XComplete = X(all(~isnan(X),2),:);
[~,C,S] = normalize(XComplete);
XNormed = normalize(X,'center',C,'scale',S);
I = isnan(X);
imputedData = knnimpute(XNormed,floor(size(XComplete,1)/2));
[U,S,V] = svd(imputedData,'econ','vector');
remained_num = size(X,2) - rank;
sigma2 = sum(S(rank+1:end).^2)./remained_num;
RegS = diag(S(1:rank)) - sigma2*diag(1./S(1:rank));
FA_Approx = U(:,1:rank)*(RegS)*V(:,1:rank)';
NewImputed = imputedData;
NewImputed(I) = FA_Approx(I);
error = norm(NewImputed(I) - imputedData(I));
while error > 1e-8
    imputedData = NewImputed;
    [U,S,V] = svd(imputedData,'econ','vector');
    remained_num = size(X,2) - rank;
    sigma2 = sum(S(rank+1:end).^2)./remained_num;
    RegS = diag(S(1:rank)) - sigma2*diag(1./S(1:rank));
    PCA_score = U(:,1:rank)*(RegS);
    FA_Approx = PCA_score*V(:,1:rank)';
    NewImputed = imputedData;
    NewImputed(I) = FA_Approx(I);
    error = norm(NewImputed(I) - imputedData(I));
end
lowApprox = NewImputed;
