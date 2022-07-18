function result = CalcMeasures(Y, predY)
% return: result = [AC, nmi_value];
if size(Y,2) ~= 1
    Y = Y';
end;
if size(predY,2) ~= 1
    predY = predY';
end;

PredLabel = predY;
% bestMap
predY = bestMap(Y, predY);
if size(Y)~=size(predY)
    predY=predY';
end
if size(Y)~=size(PredLabel)
    PredLabel=PredLabel';
end

error_cnt = sum(Y ~= predY);
AC = length(find(Y == predY))/length(Y);
[~,nmi_value,~] = compute_nmi(Y', PredLabel');

result = [AC, nmi_value, error_cnt];
