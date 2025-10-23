function Measures = GeneralNetworkMeasures(ConnectomeMatrices)
% Clustering Coefficients normalised

Nodenum = size(ConnectomeMatrices,1);
K=sum(ConnectomeMatrices~=0,2);
density = sum(K)./(Nodenum*(Nodenum-1));
W_cuberoot = ConnectomeMatrices.^(1/3);
Z = pagemtimes(pagemtimes(W_cuberoot,W_cuberoot),W_cuberoot);
diagZ = Z.*eye(Nodenum);
C_g = squeeze(sum(pagemtimes(diagZ,ones(Nodenum,1)))./sum(K.*(K-1))./density);
clear W_cuberoot Z diagZ 

%Associativity
M = sum(K);
Mask = ConnectomeMatrices>0;
j_i = sum(ConnectomeMatrices,1);
k_i = sum(ConnectomeMatrices,2);
A_numer = sum(j_i.*k_i.*Mask,[1 2])./M - (sum(0.5*(j_i+k_i).*Mask,[1 2])./M).^2;
A_denom = sum(0.5*(j_i.^2+k_i.^2).*Mask,[1 2])./M - (sum(0.5*(j_i+k_i).*Mask,[1 2])./M).^2;
A_g = squeeze(A_numer./A_denom);
clear Mask A_numer A_denom

%Strength Variance

Strength_var = squeeze(var(sum(ConnectomeMatrices)));


%Routing Eff
G_Eff = zeros(size(ConnectomeMatrices,3),1);
for parti = 1:size(ConnectomeMatrices,3)
    G_Eff(parti) = rout_efficiency(ConnectomeMatrices(:,:,parti),'inv');
end


%Diffusion Eff

Diff_Eff = zeros(size(ConnectomeMatrices,3),1);
for parti = 1:size(ConnectomeMatrices,3)
    Diff_Eff(parti) = diffusion_efficiency(ConnectomeMatrices(:,:,parti));
end


Measures = [C_g A_g Strength_var G_Eff Diff_Eff];


