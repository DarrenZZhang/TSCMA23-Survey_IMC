function Clu_result = OPIMC(X, W, option, block_size)
% If you use the code, please cite the following papers:
% [1] Hu M, Chen S. One-pass incomplete multi-view clustering[C]//Proceedings of the AAAI conference on artificial intelligence. 2019, 33(01): 3838-3845.
% [2] Jie Wen, Zheng Zhang, Lunke Fei, Bob Zhang, Yong Xu, Zhao Zhang, Jinxing Li, A Survey on Incomplete Multi-view Clustering, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS: SYSTEMS, 2022.
% Thanks Menglei Hu for providing the source code of OPIMC!


    rng('default')
    rand('seed',1817)
    num_passes = option.pass;
    num_views = numel(X);
    total = size(X{1}, 2);                  
    label = option.label;
    alpha = option.alpha;
    
    skip_loss = option.loss;
    maxIter = option.maxiter;
    
    k = option.k;                           
    tol = option.tol;
    num_feature = zeros(num_views,1);
    U = cell(num_views,1);
    for i = 1:num_views
        num_feature(i) = size(X{i},1);      
        U{i} = rand(num_feature(i), k);     
    end
    R = cell(num_views, 1);
    T = cell(num_views, 1);
    
    num_block = ceil(total/block_size);     
    Loss = zeros(num_passes, num_block);
    
    average_nmi = zeros(num_passes, 1);
    average_ac  = zeros(num_passes, 1);
    label_total = zeros(num_passes,size(X{1},2));
    for pass = 1:num_passes
        t1 = tic;   
        if pass == 1
            sum_num = 0; 
            for i = 1:num_views
                R{i} = zeros(size(U{i},1),k); 
                T{i} = zeros(k);  
            end
        else
            label_total(pass,:) = label_total(pass-1,:);   
        end
        for block_index = 1:num_block
            data_range  = (block_index-1)*block_size+1:block_index*block_size;
            V = rand(block_size,k);             
            if block_index==num_block
                data_range = (block_index-1)*block_size+1:total;
                V = rand(total-(num_block-1)*block_size,k);
            end
            % Get data blocks, W blocks.
            X_block = cell(num_views,1);
            W_block = cell(num_views,1);
            for i = 1:num_views
                X{i}(:,data_range) = NormalizeFea(X{i}(:,data_range),0);    % ·¶Êý¹éÒ»»¯
                X_block{i} = X{i}(:,data_range);
                W_block{i} = diag(W{i}(data_range));
            end
            % Get initial U and V;  fix U and update V
            if pass == 1
                sum_num = sum_num + size(X_block{1},2);
                if(block_index==1)
                    % alternating initialize U and V
                    for j = 1:num_views
                         U{j} = rand(size(U{j}));
                    end
                    V = rand(size(V));    
                    for i= 1:size(V,1)
                        V(i,:) = V(i,:)/sum(V(i,:));
                    end 
                    [~,label_total(pass,data_range)] = max(V');
                else
                    % Only initialize V
                    V = UpdateV(X_block,W_block,U,V,num_views,option);
                    [~,label_total(pass,data_range)] = max(V');
                end 
            else
                % Only use previous V
                V = full(sparse(1:length(data_range),label_total(pass,data_range),1,length(data_range),k,length(data_range)));
            end 
            iter = 0;
            converge =0;
            
            log_out = 0;
            for i = 1:num_views
                tmp1 = R{i}+ X_block{i}*V; 
                tmp2 = T{i}+ V'*W_block{i}*V;
                log_out = log_out - 2*trace(U{i}'*tmp1) + trace(U{i}'*U{i}*tmp2) + alpha*norm(U{i},'fro')^2;
            end
                            
            % update U and V 
            while(iter<maxIter && converge == 0)      
                if pass ~= 1 && iter ==0
                    V_pre = full(sparse(1:length(data_range),label_total(pass-1,data_range),1,length(data_range),k,length(data_range)));
                    for i = 1:num_views
                        T{i} = T{i} - V_pre'*W_block{i}*V_pre;      
                        R{i} = R{i} - X_block{i}*V_pre;
                    end
                end
               
                for i = 1:num_views
                    tmp1 = T{i} + V'*W_block{i}*V + alpha*eye(k);
                    tmp2 = R{i} + X_block{i}*V;                    
                    U_new = tmp2/tmp1;                            
                    if pass == 1 & block_index == 1 
                        U_new(:,find(diag(V'*W_block{i}*V)==0)) = repmat(mean(X_block{i},2),1,length(find(diag(V'*W_block{i}*V)==0)));
                        U{i} = U_new;
                    else
                        U_new(:,find(diag(T{i} + V'*W_block{i}*V)==0)) = U{i}(:,(find(diag(T{i} + V'*W_block{i}*V)==0)));
                        U{i} = U_new;
                    end
                end

                V = UpdateV(X_block,W_block,U,V,num_views,option);               
                [~,label_total(pass,data_range)] = max(V');
                 
                log_out_new = 0;
                for i = 1:num_views
                    tmp1 = R{i}+ X_block{i}*V; 
                    tmp2 = T{i}+ V'*W_block{i}*V;
                    log_out_new = log_out_new - 2*trace(U{i}'*tmp1) + trace(U{i}'*U{i}*tmp2) + alpha*norm(U{i},'fro')^2;
                end
                
                if abs((log_out_new - log_out)/log_out) < tol
                    converge = 1;
                end
                log_out = log_out_new;
                iter = iter + 1;
            end
            for i = 1:num_views
                R{i} = R{i} + X_block{i}*V;
                T{i} = T{i} + V'*W_block{i}*V;
            end
            
            if skip_loss == 0
                Loss(pass, block_index) = log_out_new/sum_num;  
            end
        end
        ttime = toc(t1); 
        result_LatLRR = ClusteringMeasure(label, label_total(pass,:)')*100;
        ac = result_LatLRR(1);
        nmi_value = result_LatLRR(2);
        purity    = result_LatLRR(3);
        average_nmi(pass, 1) = nmi_value;
        average_ac(pass, 1) = ac;
        average_pur(pass, 1) = purity;
        [ARi(pass, 1),~,~,~] = RandIndex(label, label_total(pass,:)');   
    end
    Clu_result.ACC    = average_ac(pass, 1);
    Clu_result.NMI    = average_nmi(pass, 1);
    Clu_result.Purity = average_pur(pass, 1);
    Clu_result.ARi    = ARi(pass,1);     
end
    

    