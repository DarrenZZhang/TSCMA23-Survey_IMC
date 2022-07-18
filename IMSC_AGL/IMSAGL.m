function [U] = IMSAGL2(X,G,Z_ini,F_ini,invXX,num_clust,lambda1,lambda2,lambda3,miu,rho,max_iter)
% % % The code is written by Jie Wen, if you have any problems, 
% % % please don't hesitate to contact me: jiewen_pr@126.com
% % % If you find the code is useful, please cite the following reference:
% % % Jie Wen, Yong Xu, Hong Liu, Incomplete Multi-view Spectral Clustering with
% % % Adaptive Graph Learning, IEEE Transactions on Cybenetics, 2019.
num_view = length(X);
Z = Z_ini;
S = Z;
W = Z;
F = F_ini;
for i = 1:num_view
    E{i}  = zeros(size(X{i}));
    C1{i} = zeros(size(X{i}));
    C2{i} = zeros(size(Z{i}));
    C3{i} = zeros(size(Z{i}));
end
for iter = 1:max_iter
%     fprintf('--------------iter number------------')
%     iter
    XXv = 0;
    % ------------ U -------------- %
    FFF = 0;
    for i = 1:num_view
        FFF = FFF+F{i}*F{i}';
    end
    FFF(isnan(FFF))=0;
    FFF(isinf(FFF))=1e5;
    try
        [V,D] = eig(FFF);
    catch ME
        if (strcmpi(ME.identifier,'MATLAB:eig:NoConvergence'))
            [V,D] = eig(FFF, eye(size(FFF)));
        else
            rethrow(ME);
        end
    end  
    [D_sort, ind] = sort(diag(D),'descend');
    U = V(:, ind(1:num_clust));   
%     try
%         [U,D] = eigs(FFF,num_clust,'la');
%     catch ME
%         if (strcmpi(ME.identifier,'MATLAB:eig:NoConvergence'))
%             [U,D] = eigs(FFF, eye(size(FFF)),num_clust,'la');
%         else
%             rethrow(ME);
%         end
%     end  
       
    FFF = 0; 
    linshi_obj1 = 0;
    linshi_obj2 = 0;    
    for i = 1:num_view
        % --------------1 Z{i} ------------ % 
        G1 = X{i}-E{i}+C1{i}/miu;
        G2 = S{i} - C2{i}/miu;
        G3 = W{i} - C3{i}/miu;
        Z{i} = invXX{i}*(X{i}'*G1+G2+G3);
        clear G1 G2 G3
        % ------------2 W{i} -------------- %
        P = G{i}*F{i};
        Q = L2_distance_1(P',P');
        M = Z{i}+C3{i}/miu;
        linshi_W = M-0.5*lambda1*Q/miu;
        linshi_W = linshi_W-diag(diag(linshi_W));
        for ic = 1:size(Z{i},2)
            ind = 1:size(Z{i},2);
            ind(ic) = [];
            W{i}(ic,ind) = EProjSimplex_new(linshi_W(ic,ind));
        end
        clear linshi_W P Q M ind ic
        % ----------------3 S{i} ---------------- %
        temp = Z{i}+C2{i}/miu;
        [AU,SU,VU] = svd(temp,'econ');
        AU(isnan(AU)) = 0;
        VU(isnan(VU)) = 0;
        SU(isnan(SU)) = 0;
        SU = diag(SU);    
        SVP = length(find(SU>1/miu));
        if SVP >= 1
            SU = SU(1:SVP)-1/miu;
        else
            SVP = 1;
            SU = 0;
        end
        S{i} = AU(:,1:SVP)*diag(SU)*VU(:,1:SVP)';
        clear temp AU VU SVP
        % -------------- 4 F{i} --------- %
        WW = (abs(W{i})+abs(W{i}'))*0.5;
        LL = diag(sum(WW))-WW;
        M = lambda3*U*U'-lambda1*G{i}'*LL*G{i};
        M(isnan(M))=0;
        M(isinf(M))=1e5;
        try
            [U,D] = eigs(M,num_clust,'lr');
        catch ME
            if (strcmpi(ME.identifier,'MATLAB:eig:NoConvergence'))
                [U,D] = eigs(M, eye(size(M)),num_clust,'lr');
            else
                rethrow(ME);
            end
        end  
%         [D_sort, ind] = sort(diag(D),'descend');
        F{i} = U;   
        clear M WW 
        FFF = FFF+F{i}*F{i}';
        % -------------5 E{i} -------- %
        temp1 = X{i}-X{i}*Z{i}+C1{i}/miu;
        temp2 = lambda2/miu;
        E{i} = max(0,temp1-temp2)+min(0,temp1+temp2);
        clear temp1 temp2
        % ----------- C1 C2 C3 --------- %
        leq1 = X{i}-X{i}*Z{i}-E{i};
        leq2 = Z{i}-S{i};
        leq3 = Z{i}-W{i};
        C1{i} = C1{i} + miu*leq1;
        C2{i} = C2{i} + miu*leq2;
        C3{i} = C3{i} + miu*leq3;
        % ---------- obj --------- %
        linshi_obj1 = linshi_obj1+sum(abs(SU))+lambda1*trace(F{i}'*G{i}'*LL*G{i}*F{i})+lambda2*sum(abs(E{i}(:)));
        linshi_obj2 = linshi_obj2+norm(leq1,'fro')^2+norm(leq2,'fro')^2+norm(leq3,'fro')^2;
        XXv = XXv + norm(X{i},'fro')^2;
    end
    % ----------------- obj ----------- %
    obj(iter) = (linshi_obj1+linshi_obj2+lambda3*(num_view*num_clust-trace(U'*FFF*U)))/XXv;
    clear FFF
    % ---------- miu ------------- %
    miu = min(rho*miu,1e8);
    if iter > 2 && abs(obj(iter)-obj(iter-1))<1e-7
        iter
        break;
    end
end

end