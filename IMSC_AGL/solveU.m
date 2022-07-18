function U = solveU(F,numOfClasses)

numOfViews = length(F);
M = 0;
for i = 1:numOfViews
    M = M+F{i}*F{i}';
end

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
[D_sort, ind] = sort(diag(D),'descend');
U = V(:, ind(1:numOfClasses));
end




