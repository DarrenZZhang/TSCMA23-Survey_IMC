function [Clu_result] = BSV(X,ind_folds,ground_label,numClust)
% If you use the code, please cite the following references:
% [1] Jie Wen, Zheng Zhang, Lunke Fei, Bob Zhang, Yong Xu, Zhao Zhang, Jinxing Li, A Survey on Incomplete Multi-view Clustering, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS: SYSTEMS, 2022.
% [2] Wen J, Xu Y, Liu H. Incomplete multiview spectral clustering with adaptive graph learning[J]. IEEE Transactions on Cybernetics, 2020, 50(4): 1418-1429.
% The following released codes are written  by Jie Wen. For any problems, contact jiewen_pr@126.com
 
if size(X{1},2)~=length(ground_label)
    for iv = 1:length(X)
       X{iv} = X{iv}'; 
    end
end

for iv = 1:length(X)
    X1 = X{iv};
    X1 = NormalizeFea(X1,0);
    ind_1 = find(ind_folds(:,iv) == 1);
    linshi_X = X1(:,ind_1);
    mean_linshi_X = mean(linshi_X,2);       % 一个列

    ind_0 = find(ind_folds(:,iv) == 0);
    X1(:,ind_0) = repmat(mean_linshi_X,1,length(ind_0));   %缺失视角幅值均值
    rand('seed',3523)
    pre_labels    = kmeans(X1',numClust);
    result_LatLRR = ClusteringMeasure(ground_label, pre_labels);       
    AC(iv)    = result_LatLRR(1)*100;
    NMI(iv) = result_LatLRR(2)*100;
    Purity(iv)= result_LatLRR(3)*100;
    [ARi(iv),~,~,~] = RandIndex(ground_label, pre_labels);               
end

[max_ACC,id] = max(AC);
Clu_result.ACC = max_ACC;
Clu_result.NMI = NMI(id);
Clu_result.Purity = Purity(id);
Clu_result.ARi = ARi(id);
end