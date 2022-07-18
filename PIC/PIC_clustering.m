function [Clu_result] = PIC_clustering(Dataname,percentDel,f)
% If you use the code, please cite the following papers:
% Wang H, Zong L, Liu B, et al. Spectral perturbation meets incomplete multi-view data[C]//Proceedings of the 28th International Joint Conference on Artificial Intelligence. 2019: 3677-3683.

rand('seed',5555)
opts = [];
opts.repeat = 2;
opts.normalized = 1;
opts.beta = 0.1; 
Datafold = [Dataname,'_percentDel_',num2str(percentDel),'.mat'];
load(Dataname);
load(Datafold);
num_view = length(X);
numClust = length(unique(truth));
numInst  = length(truth); 
ind_folds = folds{f};

if size(X{1},1)~=numInst
    for iv = 1:num_view
        X{iv} = X{iv}';
    end
end
view_num = length(X);
nCluster = length(unique(truth));
for iv = 1:view_num
    ind_0 = find(ind_folds(:,iv) == 0);
    X{iv}(ind_0,:) = NaN;
    ind_1 = find(ind_folds(:,iv) == 1);
    index{iv} = ind_1;
    X{iv}=X{iv}';
end
data = X;
clear X
% === Normalization =======
n_num = size(data{1},2); % the number of instances
if opts.normalized
    for i = 1:view_num
        for  j = 1:n_num
            normItem = std(data{i}(:,j));
            if (0 == normItem)
                normItem = eps;
            end;
            data{i}(:,j) = (data{i}(:,j)-mean(data{i}(:,j)))/(normItem);
        end;
    end;
end
% === Similarity matrix ========
pn = 13;
W = cell(1, view_num);
distX = cell(1, view_num);
for i = 1 :view_num
    [W{i}, distX{i}] = constructS_PNG(data{i}, pn, 0);
end
% === Similarity matrix completion ===
S = SCompletion(W, index, 1);
D = cell(1, view_num);
for v = 1:view_num
    D{v} = diag(sum(S{v},2));
end
fprintf('creating similarity matrix is done\n');

[predictLabel, weight, theta] = PIC(S, D, nCluster, view_num, opts);
result_CLU = ClusteringMeasure(truth, predictLabel)*100;   

Clu_result.ACC = result_CLU(1);
Clu_result.NMI = result_CLU(2);
Clu_result.Purity = result_CLU(3);
[Clu_result.ARi,~,~,~] = RandIndex(truth, predictLabel);


end