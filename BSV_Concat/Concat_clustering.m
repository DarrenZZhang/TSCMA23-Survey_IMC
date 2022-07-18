function [Clu_result] = Concat_clustering(Dataname,percentDel,f)
% If you use the code, please cite the following references:
% [1] Jie Wen, Zheng Zhang, Lunke Fei, Bob Zhang, Yong Xu, Zhao Zhang, Jinxing Li, A Survey on Incomplete Multi-view Clustering, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS: SYSTEMS, 2022.
% [2] Wen J, Xu Y, Liu H. Incomplete multiview spectral clustering with adaptive graph learning[J]. IEEE Transactions on Cybernetics, 2020, 50(4): 1418-1429.
% The following released codes are written  by Jie Wen. For any problems, contact jiewen_pr@126.com
 
Datafold = [Dataname,'_percentDel_',num2str(percentDel),'.mat'];
load(Dataname);
load(Datafold);
[numFold,numInst] = size(folds);
numClust = length(unique(truth));
ind_folds = folds{f};

[Clu_result] = Concat(X,ind_folds,truth,numClust);
end