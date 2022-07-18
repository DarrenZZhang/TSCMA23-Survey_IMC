function Clu_result = DAIMC_clustering(Dataname,percentDel,f)
% If you use the code, please cite the following papers:
% [1] Hu M, Chen S. Doubly aligned incomplete multi-view clustering[C]//Proceedings of the 27th International Joint Conference on Artificial Intelligence. 2018: 2262-2268.
% [2] Jie Wen, Zheng Zhang, Lunke Fei, Bob Zhang, Yong Xu, Zhao Zhang, Jinxing Li, A Survey on Incomplete Multi-view Clustering, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS: SYSTEMS, 2022.
% Thanks Menglei Hu for providing the source code of DAIMC!

% rand('seed',7878)
% rand('seed',6812)
rand('seed',6821)
Datafold = [Dataname,'_percentDel_',num2str(percentDel),'.mat'];
load(Dataname);
load(Datafold);
num_view = length(X);
numClust = length(unique(truth));
numInst  = length(truth); 
ind_folds = folds{f};

options.afa = 0.0001;
options.beta = 10000; 
% options.afa = 0.0001;
% options.beta = 100; 
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

clear X X1 W1 ind_0
X = Y;
clear Y
[U0,V0,B0] = newinit(X,W,numClust,num_view);
[U,V,B,F,P,N] = DAIMC(X,W,U0,V0,B0,truth,numClust,num_view,options);

% indic = litekmeans(V, numClust, 'Replicates', 20);
indic = kmeans(V, numClust);
result_CLU = ClusteringMeasure(truth, indic)*100;   

Clu_result.ACC = result_CLU(1);
Clu_result.NMI = result_CLU(2);
Clu_result.Purity = result_CLU(3);
[Clu_result.ARi,~,~,~] = RandIndex(truth, indic);

                
end