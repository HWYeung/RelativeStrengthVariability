function [rsquaresContribute, ConditionNumber, ...
    lnLike_before, lnLike_after, ...
    Significance] = GetR2Contributions(Covar,gfdiffusion,bagsofMeasure,...
    g,numInterest, Index)
%bagsofMeasure, a struct variable with fields: name, measures
colours = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
LBag = length(bagsofMeasure);
rsquarelengths = 2+LBag - 1;
rsquaresContribute = zeros(size(gfdiffusion,2),rsquarelengths,rsquarelengths+1);
ConditionNumber = zeros(size(gfdiffusion,2),rsquarelengths);
lnLike_before = zeros(size(gfdiffusion,2),rsquarelengths);
lnLike_after = zeros(size(gfdiffusion,2),rsquarelengths);
%stats = zeros(4,4,6,6);
if ~isempty(Index)
    gfdiffusion = gfdiffusion(Index,:);
    for i = 1:length(bagsofMeasure)
        bagsofMeasure(i).measure = bagsofMeasure(i).measure(Index,:,:);
    end
end



OneorSix = arrayfun(@(x) size(x.measure,3),bagsofMeasure); %should be 1 or 6, matching number of weights
I_1_6 = OneorSix == 6;
for i = 1:size(gfdiffusion,2)
    Indexing = max(1,I_1_6*i);
    MeasureOfInterest = cell(numInterest,1);
    L_int = numInterest - 1;
    for interest = L_int:-1:0
        MeasureOfInterest{end-interest} = bagsofMeasure(end-interest).measure(:,:,Indexing(end-interest));
    end
    X = cell(rsquarelengths-L_int,1);
    X{1} = Covar; X{2} = gfdiffusion(:,i);
    for measureset = 3:rsquarelengths-L_int
        X{measureset} = bagsofMeasure(measureset-2).measure(:,:,Indexing(measureset-2));
    end

    for modelorder = 1: rsquarelengths
        X_numbering = max([modelorder-L_int,1]);
        MOI_numbering = min(modelorder, numInterest);
        [Contributions, condition,...
            lnlike_before, lnlike_after] = R2ContriLM(X(1:X_numbering),MeasureOfInterest(1:MOI_numbering),g);
        ConditionNumber(i,modelorder) = condition;
        lnLike_before(i, modelorder) = lnlike_before;
        lnLike_after(i, modelorder) = lnlike_after;
        rsquaresContribute(i,modelorder,1:X_numbering) = Contributions(1:X_numbering);
        rsquaresContribute(i,modelorder,end-L_int:end-numInterest+MOI_numbering) = Contributions(end-MOI_numbering+1:end);
    end
end

LikelihoodDiff = lnLike_after - lnLike_before;
Significance = arrayfun(@(x) chi2cdf(x,1,'upper'),LikelihoodDiff);

R2stacks = cumsum(rsquaresContribute,3);
%rsquaresContribute = R2stacks;
labels = cell(rsquarelengths+1,1);
labels{1} = "Covariates"; labels{2} = "Mean Edge Weight";
for measureset = 3:rsquarelengths+1
    labels{measureset}=bagsofMeasure(measureset-2).name;
end

minimums = min(R2stacks,[],'all') - 0.0005;
maximums = max(R2stacks,[],'all') + 0.0005;
%R2stacks = R2stacks(:,:,3:end);
%colours = {[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],[0 0.4470 0.7410]};
for stacking = 1:rsquarelengths+1
    bar(1:6,R2stacks(:,:,rsquarelengths+2-stacking),'FaceColor',colours{rsquarelengths+2-stacking});
    hold on
end
set(gca,'FontSize',12)
ylabel('Rsquared Contributions')
xlabel('Connectome Weights')
set(gca,'XTickLabel',{'MD','FA','SC','OD','ISOVF','ICVF'});
title('Model fit for predicting the cognitive g-factor')
hBLG = bar(nan(2,rsquarelengths+1));         % the bar object array for legend
hold on
for i=1:rsquarelengths+1
  hBLG(i).FaceColor=colours{i};
end
hLG=legend(hBLG,labels,'location','eastoutside');
hold off
ylim([minimums maximums])

figure
R2stacks = cumsum(rsquaresContribute(:,:,2:end),3);
labels(1) = [];
minimums = min(R2stacks,[],'all');
maximums = max(R2stacks,[],'all') + 0.0005;
for stacking = 1:rsquarelengths
    bar(1:6,R2stacks(:,:,rsquarelengths+1-stacking),'FaceColor',colours{rsquarelengths+1-stacking});
    hold on
end
set(gca,'FontSize',12)
ylabel('Incremental Rsquared Contributions')
xlabel('Connectome Weights')
set(gca,'XTickLabel',{'MD','FA','SC','OD','ISOVF','ICVF'});
title('Model fit for predicting the cognitive g-factor')
hBLG = bar(nan(2,rsquarelengths));         % the bar object array for legend
hold on
for i=1:rsquarelengths
  hBLG(i).FaceColor=colours{i};
end
hLG=legend(hBLG,labels,'location','eastoutside');
hold off
ylim([minimums maximums])



end

function [Contributions, cond_num, lnlike_before, ...
    lnlike_after] = R2ContriLM(X,Interest,y)
Xmat = cat(2,X{:});
Groups = [];
for mea_size = 1:length(X)
    Groups = [Groups ; mea_size*ones(size(X{mea_size},2),1)];
end
InterestMat = cat(2,Interest{:});
for MOI_size = 1:length(Interest)
    Groups = [Groups ; (MOI_size+length(X))*ones(size(Interest{MOI_size},2),1)];
end
[~,R] = qr(normalize([Xmat InterestMat]),0);
cond_num = rcond(R);
[X_resid] = GSRegressResid(normalize([Xmat InterestMat]));
correlations = corr(X_resid,y);
mdl = fitglm(X_resid,y);
b = mdl.Coefficients.Estimate(2:end);
temptContri = b.*correlations;
Contributions = splitapply(@sum,temptContri,Groups);
lnlike_after = mdl.LogLikelihood;
mdl = fitglm(X_resid(:,1:end-1),y);
lnlike_before = mdl.LogLikelihood;
end

