clc
clear all
close all
%% Transition matrix for hierarchical clustering

selectTechnq = 1;                                                          % 1 for GNG , 2 for Kmeans

if selectTechnq == 1
    load('VocabGNG.mat')
elseif selectTechnq == 2
    load('VocabKMeans.mat')
elseif selectTechnq == 3
    load('VocabDBSCAN.mat')
else
    load('VocabSOM.mat')
    
end
% Load hierarchical clustering
load('VocabHier.mat')
% Load data
load('PMDatafile.mat')

InputData =[structSyncData.Filtered.xPos, structSyncData.Filtered.yPos,...
    structSyncData.Filtered.divxPos, structSyncData.Filtered.divyPos ];

cluster =[];

%  Transition matrix
for i = 1:1:size(Hier.DataNodesLevel,2)
    tempDataNodes =  Hier.DataNodesLevel{1,i};
    maxNodes = max(Hier.labels_in_level{1, i+1});
    transitionMat = zeros(maxNodes,maxNodes);
    timesSpent = [];
    currLength = size(InputData,1) ;
    nodesInTime = Hier.ColorDataLevel{i};                                   %%%%% COMPUTE FOR EACH LEVEL
    ind = find(diff(nodesInTime) ~= 0);
    
    %% transition Matrix
    for k = 1:currLength-1
        transitionMat(nodesInTime(k,1),nodesInTime(k+1,1)) =...
            transitionMat(nodesInTime(k,1),nodesInTime(k+1,1)) + 1;
    end
    
    transitionMatCombine1 = transitionMat(any(transitionMat, 2), :);  % remove extra rows
    transitionMatCombine =  transitionMatCombine1(:, any ( transitionMatCombine1, 1));                      % remove extra columns
    
    %% time Stamp Matrix
    
    codeInd = [0; ind];
    tspentTran = diff(codeInd);
    for k = 1:size(tspentTran,1)
        if size(unique([timesSpent;tspentTran(k,1)]),1) ~= size(unique(timesSpent),1)
            timeMats{1,tspentTran(k)} = zeros(maxNodes,maxNodes);
            timesSpent = [timesSpent; tspentTran(k)];
        end
        timeMats{1,tspentTran(k)}(nodesInTime(ind(k),1),nodesInTime(ind(k)+1,1)) =...
            timeMats{1,tspentTran(k)}(nodesInTime(ind(k),1),nodesInTime(ind(k)+1,1)) + 1;
    end
    
    %
    
    ind2 = find(diff(nodesInTime) == 0);
    tspentSame = 1;
    for k = 1:size(ind2,1)
        if k > 1
            if ind2(k) == ind2(k-1) + 1
                tspentSame = tspentSame + 1;
            else
                tspentSame = 1;
            end
        end
        
        if size(unique([timesSpent;tspentSame]),1) ~= size(unique(timesSpent),1)
            timeMats{1,tspentSame} = zeros(maxNodes,maxNodes);
            timesSpent = [timesSpent; tspentSame];
        end
        
        timeMats{1,tspentSame}(nodesInTime(ind2(k),1),nodesInTime(ind2(k)+1,1)) =...
            timeMats{1,tspentSame}(nodesInTime(ind2(k),1),nodesInTime(ind2(k)+1,1)) + 1;
        
    end
    cluster = [cluster ; Hier.Clusters{1,i}];
    for ij=1:1:size(timeMats,2)
        
        c1 = cluster(:,1);
        c2 = cluster(:,2);
        C = [c1;c2];
        tempTimeMat = timeMats{1,ij};
        tempTimeMat(C,:)=[];
        tempTimeMat(:,C)=[];
        
        TimeMat{1,ij} = tempTimeMat;
        
    end
    
    %     cluster
    transitionMatCombine = transitionMatCombine./repmat(sum(transitionMatCombine,2) + (sum(transitionMatCombine,2)==0),1,size(transitionMatCombine,2));
    
    timeMatsLevel{1,i} = TimeMat;
    transitionMatLevel{1,i} = transitionMatCombine;
    
    
end
Hier.InputData = InputData;
Hier.transitionMatLevel = transitionMatLevel;
Hier.timeMatsLevel = timeMatsLevel;

save('VocabHier.mat','Hier')






