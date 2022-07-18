function V = UpdateV(X,W,U,V,viewNum,option)
[n,r] = size(V);
D = zeros(n,r); 
for i = 1:viewNum
    bb = full(sum(U{i}.*U{i},1));
    ab = full(X{i}'*U{i});
    D = D + W{i}*(bb(ones(1,n),:) - 2*ab);
end
[~,label] = min(D,[],2);
 V = full(sparse(1:n,label,1,n,r,n));
end


