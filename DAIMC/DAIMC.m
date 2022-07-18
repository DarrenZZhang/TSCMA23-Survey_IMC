function [U,V,B,F,P,N] = DAIMC(X,W,U,V,B,label,r,viewNum,options)
% If you use the code, please cite the following papers:
% [1] Hu M, Chen S. Doubly aligned incomplete multi-view clustering[C]//Proceedings of the 27th International Joint Conference on Artificial Intelligence. 2018: 2262-2268.
% [2] Jie Wen, Zheng Zhang, Lunke Fei, Bob Zhang, Yong Xu, Zhao Zhang, Jinxing Li, A Survey on Incomplete Multi-view Clustering, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS: SYSTEMS, 2022.
% Thanks Menglei Hu for providing the source code of DAIMC!

eta = 1e-10;
F = 0;
P = 0;
N = 0;
D = cell(viewNum,1);
for i = 1:viewNum
    for k = 1:size(B{i},1)
        D{i}(k,k) = 1/sqrt(norm(B{i}(k,:),2).^2+eta);    
    end
end

time = 0;
f =0;
while 1
    time = time+1;    
    for i = 1:viewNum
        tmp1 = options.afa*B{i}*B{i}';
        tmp2 = V'*W{i}*V;
        tmp3 = X{i}*W{i}*V + options.afa*B{i};
        U{i} = lyap(tmp1, tmp2, -tmp3);
    end  
    V = UpdateV_DAIMC(X,W,U,V,viewNum);
    Q = diag(ones(1,size(V,1))*V);
    V = V * inv(Q);
    for i = 1:viewNum
        U{i} = U{i}*Q;
        invD = diag(1./diag(0.5*options.beta*D{i}));
        B{i} = (invD - invD * U{i} * inv(U{i}'*invD*U{i} + eye(r)) * U{i}' * invD)*U{i};
       for k = 1:size(B{i},1)
           D{i}(k,k) = 1/sqrt(norm(B{i}(k,:),2).^2+eta);   
       end
    end
       
    ff = 0;
    for i = 1:viewNum
        tmp1 = (X{i} - U{i}*V')*W{i};
        tmp2 = B{i}'*U{i} - eye(r);
        tmp3 = sum(1./diag(D{i}));
        ff = ff + sum(sum(tmp1.^2)) + options.afa*(sum(sum(tmp2.^2)) + options.beta*tmp3);
    end
    F(time) = ff;  
    if abs(ff-f)/f < 1e-4 | abs(ff-f) >  1e100 | time == 30
        break;
    end
    f = ff;
end