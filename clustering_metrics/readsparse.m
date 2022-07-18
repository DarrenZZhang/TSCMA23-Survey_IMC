function [Xf] = readsparse(Xs)
%% Input: Xs: sparse format, the first line is the (numRow, numCol, numEdge). following lines are (rowid, colid, value)
%% Output: Xf, full format matrix 
 
numRow=Xs(1,1);
numCol=Xs(1,2);
Xf=zeros(numRow,numCol);
for i=2:size(Xs,1)
    Xf(Xs(i,1), Xs(i,2))=Xs(i,3);       
end