function [U,V,B] = newinit(X,W,r,viewNum)
B = cell(viewNum,1);
U = cell(viewNum,1);
H = cell(viewNum,1);
XX = cell(viewNum,1);

for i = 1:viewNum
    item = diag(W{i});
    temp = find(item == 0);
    XX{i} = X{i};
    XX{i}(:,temp)= [];
    Mx = mean(XX{i},2);
    X{i}(:,temp)= repmat(Mx,1,length(temp));
end


sumH = 0;
for i = 1:viewNum
    [d,n] = size(X{i});
    [ilabels,C] = litekmeans(X{i}', r, 'Replicates', 20);
    U{i} = C' + 0.1*ones(d,r);  
    G = zeros(n,r);
    for j=1:r
        G(:,j)=(ilabels == j*ones(n,1));
    end 
    H{i}=G+0.1*ones(n,r);
    sumH = sumH + H{i};    
end
V = sumH/viewNum;
Q = diag(ones(1,size(V,1))*V);
V = V * inv(Q);
for i = 1:viewNum
    U{i} = U{i}*Q;
end
lamda = 1e-5;
for i = 1:viewNum
    [d,~] = size(U{i});
    invI = diag(1./diag(lamda*eye(d)));
    B{i} = (invI - invI * U{i} * inv(U{i}'*invI*U{i} + eye(r)) * U{i}' * invI) * U{i};
end
end


