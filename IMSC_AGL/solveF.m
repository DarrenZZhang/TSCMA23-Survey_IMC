function F = solveF(Z,G,numOfClasses)

numOfViews = length(Z);

for i=1:numOfViews
    W = (abs(Z{i})+abs(Z{i}'))/2;    
    LapN = diag(sum(W))-W;
    M = G{i}'*LapN*G{i};
    M(isnan(M))=0;
    M(isinf(M))=1e5;
    try
        [V,D] = eig(M);
    catch ME
        if (strcmpi(ME.identifier,'MATLAB:eig:NoConvergence'))
            [V,D] = eig(M, eye(size(M)));
        else
            rethrow(ME);
        end
    end
    [D_sort, ind] = sort(diag(D));
    ind2 = find(D_sort>1e-6);
    if length(ind2)>numOfClasses
        F{i} = V(:, ind2(1:numOfClasses));
    else
        F{i} = V(:, ind(1:numOfClasses));
    end
end



