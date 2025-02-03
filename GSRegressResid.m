function U = GSRegressResid(V)
    sz = size(V);
    V = V - mean(V);
    U = zeros(sz);
    U(:,1) = V(:,1) / norm(V(:,1));
    for xpos = 2:sz(2)
        U(:,xpos) = V(:,xpos);       
        U(:,xpos) = U(:,xpos) - U(:,1:xpos-1)*(U(:,1:xpos-1)'*U(:,xpos));
        U(:,xpos) = U(:,xpos) / norm(U(:,xpos));
    end
    U = normalize(U);
end