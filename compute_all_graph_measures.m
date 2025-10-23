%Compute everything and save
disp("Extracting Connectomes")
ConnectomeExtraction

disp("Compute Node Variability and other graph measures")
Weight = ["MD", "FA", "SC", "OD", "ISOVF", "ICVF"];
weight = 1;
filename = strcat("UKB_37k_connectome_Raw_",Weight(weight),".mat");
load(filename)
filelength = length(Connectome.(Weight(weight)));
mask = readmatrix('mask_threshold_PT50.csv');
RSV = zeros(filelength,1,6);
RSV_4 = zeros(filelength,1,6);
RSV_21 = zeros(filelength,1,6);
GraphMeasures = zeros(filelength,5,6);
mean_edge_weight = zeros(filelength,6);
for weight = 1:6
    bartitle = strcat("ConnectomeWeight: ", Weight(weight)," , progress: StrengthVarability");
    disp(strcat("For ",bartitle));
    load(strcat("UKB_37k_connectome_Raw_",Weight(weight),".mat"));
    W_prior = Connectome.(Weight(weight));
    clear Connectome
    W = W_prior.*mask;
    for sample = 1:filelength
        tempt_w = W(:,:,sample);
        mean_edge_weight(sample,weight) = mean(tempt_w(tempt_w>0));
    end
    % max value normlisation
    W = W./max(A,[],[1 2]);
    clear W_prior
    [~,StrengthVar] = NodeRelStrengthVariability(W, 0, 1, 0);
    RSV(:,1,weight) = StrengthVar;
    [~, ~, windowed] = NodeRelStrengthVariability(W, 0, 1, 0, 1/4);
    RSV_21(:,1,weight) = mean(windowed, 2);
    [~, ~, windowed] = NodeRelStrengthVariability(W, 0, 1, 0, 1/20);
    RSV_4(:,1,weight) = mean(windowed, 2);
    clear W
    bartitle = strcat("ConnectomeWeight: ", Weight(weight)," , progress: Network Measures");
    disp(strcat("For ",bartitle));
    % graph measures arranged in the following order
    % clustering coeff, assortativity, strength variance, routing
    % efficiency, diffusion efficiency
    GraphMeasures(:,:,weight) = GeneralNetworkMeasures(W);
    clear NodeRel StrengthVar
end
filename = strcat("all_graph_measures_and_RSVs.mat");
save(filename, "RSV", "RSV_4", "RSV_21", "GraphMeasures", "mean_edge_weight");
