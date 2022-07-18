function [ind] = split_2view(X,Nsamp,pairPortion)
Nview = length(X);
ind = ones(Nsamp, Nview);      % 针对 一行是一个样本的数据库
paired_samples = round(pairPortion*Nsamp);
unpaired_num = round((Nsamp-paired_samples)*0.5);
rand_ind = randperm(Nsamp);
ind(rand_ind(paired_samples+(1:unpaired_num)),1) = 0;
ind(rand_ind(paired_samples+unpaired_num+1:end),2) = 0;
end