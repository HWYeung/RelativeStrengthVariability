
load("UKB_preprocessed_tabular_data.mat");
Weight = ["MD","FA","SC","OD","ISOVF","ICVF"];
load("all_graph_measures_and_RSVs.mat")
AllMeasures = zeros(length(RSV), 9, 6);
AllMeasures(:, 1, :) = RSV;
AllMeasures(:, 2, :) = RSV_4;
AllMeasures(:, 3, :) = RSV_21;
AllMeasures(:, 4, :) = mean_edge_weight;
AllMeasures(:, 5:end, :) = GraphMeasures;
Betas_Std_p = zeros(6,3,9);
All_Variables = ["RSV", "RSV_4", "RSV_21", "mean edge weight", "clustering coefficient",...
    "assortativity", "strength variance", "routing efficiency", "diffusion efficiency"];
for vars = 1:9
    for measure = 1:6
        % Non need to include site 3 as site 1, site 2 with the intercept will
        % span the space represented by site 3 one-hot encoding and cause
        % instability
        T1 = Covars(:,["Age", "Sex", "MRI_site1", "MRI_site2", "MRI_headpos_x", ...
            "MRI_headpos_y", "MRI_headpos_z"]);
        variable = squeeze(AllMeasures(CovartoCog, measure,vars));
        T6 = g_factor;
        mdl = fitglm(normalize([T1.Variables variable]),g_factor);
        Betas_Std_p(measure,1, vars) = mdl.Coefficients.Estimate(end);
        Betas_Std_p(measure,2, vars) = mdl.Coefficients.SE(end);
        Betas_Std_p(measure,3, vars) = mdl.Coefficients.pValue(end);
    end
end
[~, ~, ~, adj_p]=fdr_bh(Betas_Std_p(:,3,:));
Betas_Std_p(:,3,:) = adj_p;
subtables = {};
for vars = 1:9
    Table = array2table(squeeze(Betas_Std_p(:,:,vars)));
    Table.Properties.RowNames = Weight;
    Table.Properties.VariableNames = ["beta", "standard error", "adjusted p-values"];
    subtables{vars} = Table;
end
main_statistics_table = table(All_Variables', subtables', 'VariableNames', {'Graph Measure', 'Statistics'});
% so the p_values are page-wise
save('graph_measures_g_stats.mat', 'main_statistics_table');