
function [f1, f2, f] = OptFunctions(datanodes, InputData)

% Normalisation of data and data nodes
% find max and min of input data, for normalisation
maxData = max(InputData,[], 1);
minData = min(InputData,[], 1);

Ndata = size(InputData, 1);

dataNorm = InputData - minData;
dataNorm = dataNorm./repmat(maxData-minData, Ndata,1);

varianceData = var(dataNorm, 1);
meanVarianceDataX = (varianceData(1, 1) + varianceData(1, 2))/2;
meanVarianceDataY = (varianceData(1, 3) + varianceData(1, 4))/2;

Nclusters = size(datanodes, 2);

%Normalising the data in each cluster
datanodesNorm = cell(1, Nclusters);
for c = 1: Nclusters
    datanodesNorm{1, c} = datanodes{1, c} - minData;
    datanodesNorm{1, c} = datanodesNorm{1, c}./repmat(maxData-minData, size(datanodes{1, c}, 1),1);
end

%% Cluster variance calculation
%Finding the variance of each cluster
variancePerCluster = zeros(Nclusters, 4);
for c = 1: Nclusters 
    variancePerCluster(c, :) = var(datanodesNorm{1, c}, 0, 1);
end

for c = 1:Nclusters
f1_1(c, 1) = ReLuF1(variancePerCluster(c,1), meanVarianceDataX) +...
        ReLuF1(variancePerCluster(c,2), meanVarianceDataX);
f2_1(c, 1) = ReLuF2(variancePerCluster(c,3), meanVarianceDataY) +...
    ReLuF2(variancePerCluster(c,4), meanVarianceDataY);
end

% a high f1 is an indicator of a cluster with a high space variance
f1 = sum(f1_1)./(Nclusters*2);
% a high f2 is an indicator of a cluster with a low speed variance
f2 = sum(f2_1)./(Nclusters*2);


f =( sum(f1_1) + sum(f2_1))/(Nclusters*2);

end

