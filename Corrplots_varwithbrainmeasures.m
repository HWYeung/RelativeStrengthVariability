load("UKB_preprocessed_tabular_data.mat");
load("all_graph_measures_and_RSVs.mat");
Y_axis = "Variability Measures";
X_axis = ["Variability Measures", "Demographic Variables",...
    "Global Brain Measures", "Mean Edge Weights"];
YLabel = ["MD", "FA", "SC", "OD", "ISOVF", "ICVF"];
XLabel = {};
XLabel{1} = ["MD", "FA", "SC", "OD", "ISOVF", "ICVF"];
XLabel{2} = ["Age", "Sex"];
XLabel{3} = ["TBV", "GMV", "WMV", "Atrophy_{ICV corrected}"];
XLabel{4} = "Mean" + ["MD", "FA", "SC", "OD", "ISOVF", "ICVF"];
variability_measures = zeros(length(RSV),3,6);
variability_measures(:,1,:) = RSV;
variability_measures(:,2,:) = RSV_4;
variability_measures(:,3,:) = RSV_21;
all_correlations = cell(3,4);
all_pvals_corr = cell(3,4);
variability_name = ["RSV", "RSV_4", "RSV_21"];
for var_mea = 1:3
    f = figure('Position',[50 50 1200 1000]);
    Var = squeeze(variability_measures(:,var_mea,:));
    Xvar = {};
    Xvar{1} = Var;
    Xvar{2} = Covars(:,1:2).Variables;
    Xvar{3} = BrainVariables.Variables;
    Xvar{4} = mean_edge_weight;
    Correlations = {};
    pvals = {};
    for i = 1:4
        subtightplot(2,2,i,[0.15,0.15], [0.15 0.05],[0.1 0.1])
        XtemptVar = Xvar{i};
        YtemptVar = Var;
        YtemptLabel = YLabel;
        if size(XtemptVar,1) < size(YtemptVar,1)
            YtemptVar = YtemptVar(CovartoCog,:);
        end
        [h, Corrs, pval] = CorrHeatMap(YtemptVar,XtemptVar);
        all_correlations{var_mea, i} = Corrs;
        all_pvals_corr{var_mea, i} = p_vals;
        h.YDisplayLabels = YtemptLabel;
        h.XDisplayLabels = XLabel{i};
        h.XLabel = X_axis(i);
        h.YLabel = Y_axis;
        set(gca,"FontSize",10)
    end
    filename = sprintf("CorrelationPlot_VarMeasureWithOthers_%s.png", variability_name(var_mea));
    exportgraphics(f, filename, 'Resolution', 300);
    close(f);
end