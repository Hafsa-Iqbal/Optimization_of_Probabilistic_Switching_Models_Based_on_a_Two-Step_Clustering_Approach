clc
clear all
close all
% load data
 load('PMDatafile.mat')

InputData =[structSyncData.Filtered.xPos, structSyncData.Filtered.yPos,...
    structSyncData.Filtered.divxPos, structSyncData.Filtered.divyPos ];

meanDataNorm = mean(InputData);
stdDataNorm = std(InputData);
dataNorm = InputData - repmat(meanDataNorm,size(InputData,1),1);
inputNorm = dataNorm./stdDataNorm;


opts = statset('Display','final');
% select the number of clusters
Clusters = 23;
%% k means 
lossfunction =[];
for i =1:1:20
    
nData = size(inputNorm,1);       
rng(i)
inputNorm = inputNorm(randperm(nData), :);

opts = statset('Display','final');
[idx,C] = kmeans(inputNorm,Clusters,'Distance','cityblock','Replicates',10,'Options',opts); %

netL = vocabularydataNodes(C,idx,inputNorm,InputData);
netL.w = C;
netL.dataColorNode  = idx;

F_l = LossFunction(netL.datanodes);
lossfunction = [lossfunction ; F_l ];
dataInfo{1,i} = netL; 
i = i+1;

end

[IndxMin , minLoss] =  min(lossfunction);
[IndxMax , maxLoss] =  max(lossfunction);
 
 BestRelevant = dataInfo{1,minLoss};
 WorstRelevant = dataInfo{1,maxLoss};
 
figure;
hold on
plot(BestRelevant.w(:,1),BestRelevant.w(:,2),'kx','MarkerSize',15,'LineWidth',3) 
title 'K-MEANS'
hold off
% save data
net = BestRelevant ;
save('VocabKMeans.mat','net')
