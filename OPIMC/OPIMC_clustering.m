function Clu_result = OPIMC_clustering(Dataname,percentDel,f)
% If you use the code, please cite the following papers:
% [1] Hu M, Chen S. One-pass incomplete multi-view clustering[C]//Proceedings of the AAAI conference on artificial intelligence. 2019, 33(01): 3838-3845.
% [2] Jie Wen, Zheng Zhang, Lunke Fei, Bob Zhang, Yong Xu, Zhao Zhang, Jinxing Li, A Survey on Incomplete Multi-view Clustering, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS: SYSTEMS, 2022.
% Thanks Menglei Hu for providing the source code of OPIMC!
Datafold = [Dataname,'_percentDel_',num2str(percentDel),'.mat'];
load(Dataname);
load(Datafold);
num_view = length(X);
numClust = length(unique(truth));
numInst  = length(truth); 
ind_folds = folds{f};

if size(X{1},2)~=numInst
    for iv = 1:num_view
        X{iv} = X{iv}';
    end
end
for iv = 1:length(X)
    X1 = X{iv};
    X1 = NormalizeFea(X1,0);
    ind_0 = find(ind_folds(:,iv) == 0);
    X1(:,ind_0) = 0 ;
    Y{iv} = X1; 
    W{iv} = diag(ind_folds(:,iv));                       
end

label = truth;
ind = ind_folds;
index = randperm(length(label));
for i = 1:num_view
%     X{i} = X{i}';           
    X{i} = X{i}(:,index);   
    W{i} = ind(index,i);    
end 
label = label(index);
addpath('OPIMC')
block_size = 30;
option.label = label;
option.k       = numClust;
option.maxiter = 30;
option.tol     = 1e-6;
option.num_cluster = numClust;
option.pass = 1;
option.loss = 0;
option.alpha = 10;
[Clu_result] = OPIMC(X, W, option, block_size);
end