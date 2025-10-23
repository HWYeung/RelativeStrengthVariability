load("all_graph_measures_and_RSVs.mat")
AllMeasures = zeros(length(RSV), 8, 6);
AllMeasures(:, 1, :) = RSV;
AllMeasures(:, 2, :) = RSV_4;
AllMeasures(:, 3, :) = RSV_21;
AllMeasures(:, 4:end, :) = GraphMeasures;
permutations = [1:4 6 7 5 8];
Names = ["RSV", "RSV_4", "RSV_21", "MeanEdge", "\gamma", "r_{assort}", "V", "E_{rout}", "E_{diff}"];
Connectomes = ["(a) MD", "(b) FA","(c) SC","(d) OD","(e) ISOVF", "(f) ICVF"];
f = figure('Position',[50 50 1200 900]);
for i = 1:6
    temp_i = i;
    subtightplot(3,2,i,[0.14,0.15], [0.15 0.05],[0.1 0.1])
    [h, Corrs] = CorrHeatMap(AllMeasures(:,permutations,i));
    h.YDisplayLabels = Names(permutations);
    h.XDisplayLabels = Names(permutations);
    h.XLabel = Connectomes(i);
    set(gca,"FontSize",9)
end
filename = "CorrelationPlot_among_graph_measures.png";
exportgraphics(f, filename, 'Resolution', 300);
close(f);