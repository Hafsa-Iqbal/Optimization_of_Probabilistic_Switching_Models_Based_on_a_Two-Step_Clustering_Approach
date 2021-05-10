%
clc;
clear;
close all;

%% Load Data

load('PMDatafile.mat')
InputData =[structSyncData.Filtered.xPos, structSyncData.Filtered.yPos,...
    structSyncData.Filtered.divxPos, structSyncData.Filtered.divyPos ];

X = InputData;


% DBSCAN Clustering Algorithm
epsilon = 0.7;
MinPts = 20;
IDX = DBSCAN(X,epsilon,MinPts);

% Plot Results
PlotClusterinResult(X, IDX);
title(['DBSCAN Clustering (\epsilon = ' num2str(epsilon) ', MinPts = ' num2str(MinPts) ')']);

% centroid of clusters
dataColorNode = IDX+1;
% Export Results data samples in nodes
datanodes = cell(1,max(dataColorNode));
nData = size(InputData,1);
for c = 1:nData
    x = dataColorNode(c);
    newSelectedData= InputData(c,:);                                        %   Organize distances between nodes and the first data point in an ascending order
    datanodes{1,x} = [datanodes{1,x}; newSelectedData];
end

%  Computation of mean , covariances  and radius of generated nodes
nodesMean = [];

for i = 1:size(datanodes,2)
    %   Calculation of mean values
    nodesMean = [nodesMean; mean(datanodes{1,i})];
    %   Calculation of covariance values
    nodesCov{1,i} = cov(datanodes{1,i});
    %   Calculation of radius of acceptances
    nodesRadAccept(1,i) = sqrt(sum((3*std(datanodes{1,i})).^2));
    
end
netD.w = nodesMean;
netD.nodesMean = nodesMean;
netD.nodesCov = nodesCov;
netD.nodesRadAccept = nodesRadAccept;
netD.datanodes = datanodes;
netD.dataColorNode = dataColorNode;
netD.data = InputData;
save('VocabDBSCAN.mat','netD')