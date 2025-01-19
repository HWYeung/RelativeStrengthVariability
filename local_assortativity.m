function [Local_assort] = local_assortativity(W)
% W = R x R x N matrix where R is the number of nodes in Adjacency matrix and 
%                            N is the number of participants

M = sum(W,[1 2]);
Mask = W;
j_i = sum(W,1);
k_i = sum(W,2);
fixed_top = (sum(0.5*(j_i+k_i).*Mask,[1 2])./M).^2;
fixed_bottom = sum(0.5*(j_i.^2+k_i.^2).*Mask,[1 2])./M - (sum(0.5*(j_i+k_i).*Mask,[1 2])./M).^2;
Varyingpart = j_i.*k_i.*Mask;
numer = (sum(Varyingpart)-j_i.*fixed_top)./M;
Local_assort = squeeze(numer./fixed_bottom);
