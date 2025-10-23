%Grouping for contribution plot

load("UKB_preprocessed_tabular_data.mat");
load("all_graph_measures_and_RSVs.mat");
bag_of_measures = struct();
bag_of_measures(1).name = "Network Measures";
bag_of_measures(1).measure = GraphMeasures;
bag_of_measures(2).name = "variability_{weights}";
bag_of_measures(2).measure = RSV;
bag_of_measures(3).name = "variability_{sliding window 21}";
bag_of_measures(3).measure = RSV_21;
bag_of_measures(4).name = "variability_{sliding window 4}";
bag_of_measures(4).measure = RSV_4;

f = figure('Position',[50 50 1500 1000]);
[rsquaresContribute, ConditionNumber, ...\
    lnLike_before, lnLike_after, ...
    Significance] = GetR2Contributions(Covars.Variables, ...
    mean_edge_weight,bag_of_measures([1 2]),g_factor,1,CovartoCog);
filename = "incremental_variance_rsv.png";
exportgraphics(f, filename, 'Resolution', 300);
close(f);

f = figure('Position',[50 50 1500 1000]);
[rsquaresContribute, ConditionNumber, ...\
    lnLike_before, lnLike_after, ...
    Significance] = GetR2Contributions(Covars.Variables, ...
    mean_edge_weight,bag_of_measures([1:3]),g_factor,2,CovartoCog);
filename = "incremental_variance_rsv_slide_21.png";
exportgraphics(f, filename, 'Resolution', 300);
close(f);

f = figure('Position',[50 50 1500 1000]);
[rsquaresContribute, ConditionNumber, ...\
    lnLike_before, lnLike_after, ...
    Significance] = GetR2Contributions(Covars.Variables, ...
    mean_edge_weight,bag_of_measures([1 2 4]),g_factor,2,CovartoCog);
filename = "incremental_variance_rsv_slide_4.png";
exportgraphics(f, filename, 'Resolution', 300);
close(f);