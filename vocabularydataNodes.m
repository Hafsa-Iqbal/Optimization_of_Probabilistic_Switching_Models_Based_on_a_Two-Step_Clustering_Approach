function netL = vocabularydataNodes(centroids,idx,inputNorm,InputData)

%Export data samples in nodes
datanodes = cell(1,max(idx));
datanodesNorm = cell(1,max(idx));
dataColorNode = [];
nData = size(inputNorm,1);
for c = 1:nData       
    
    x = inputNorm(c,:);
    d = pdist2(x, centroids,'euclidean');                                          %   pdist- Pairwise distance between two sets of observations(Eucledian distance between input and 2 nodes initialised before)
    [~, minNode] = min(d);                                                 %   Organize distances between nodes and the first data point in an ascending order
     dataColorNode = [dataColorNode; minNode];
     
     datanodesNorm{1,minNode} = [datanodesNorm{1,minNode}; x];
     
     xd = InputData(c,:);
     
    datanodes{1,minNode} = [datanodes{1,minNode}; xd];                      % normalize and ordered data

end

%%
nodesMean = [];
nodesMeanNorm = [];
for i = 1:size(datanodes,2)
    %   Calculation of mean values
    nodesMean = [nodesMean; mean(datanodes{1,i})];
    nodesMeanNorm = [nodesMeanNorm; mean(datanodesNorm{1,i})];
    
    %   Calculation of covariance values
    nodesCov{1,i} = cov(datanodes{1,i});
    nodesCovNorm{1,i} = cov(datanodesNorm{1,i});
    
    %   Calculation of radius of acceptances
    tempData = datanodes{1,i};
    nodesRadAccept_State(1,i) = sqrt(sum((3*std(tempData(:,1:2))).^2));
    nodesRadAccept_Div(1,i) = sqrt(sum((3*std(tempData(:,3:4))).^2));
    nodesRadAcceptNorm(1,i) =sqrt(sum((3*std(datanodesNorm{1,i})).^2));

 
    %   Calculation of radius of acceptances
    nodesRadAccept(1,i) = sqrt(sum((3*std(datanodes{1,i})).^2));
    
end

%%  Transition matrix

transitionMat = zeros(size(centroids,1),size(centroids,1));
timesSpent = [];
currLength = size(InputData,1) ;
nodesInTime = dataColorNode;
ind = find(diff(nodesInTime) ~= 0);

%% transition Matrix
for k = 1:currLength-1
    transitionMat(nodesInTime(k,1),nodesInTime(k+1,1)) =...
        transitionMat(nodesInTime(k,1),nodesInTime(k+1,1)) + 1;
end

%% time Stamp Matrix

codeInd = [0; ind];
tspentTran = diff(codeInd);

for k = 1:size(tspentTran,1)
    if size(unique([timesSpent;tspentTran(k,1)]),1) ~= size(unique(timesSpent),1)
        timeMats{1,tspentTran(k)} = zeros(size(centroids,1),size(centroids,1));
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
        timeMats{1,tspentSame} = zeros(size(centroids,1),size(centroids,1));
        timesSpent = [timesSpent; tspentSame];
    end
    
    timeMats{1,tspentSame}(nodesInTime(ind2(k),1),nodesInTime(ind2(k)+1,1)) =...
        timeMats{1,tspentSame}(nodesInTime(ind2(k),1),nodesInTime(ind2(k)+1,1)) + 1;
end

transitionMat = transitionMat./repmat(sum(transitionMat,2) + (sum(transitionMat,2)==0),1,size(centroids,1));

netL.transitionMat = transitionMat;
netL.TimeMats = timeMats;




%% 
netL.wNorm = centroids;
netL.nodesMean = nodesMean;
netL.nodesMeanNorm = nodesMeanNorm;
netL.nodesCov = nodesCov;
netL.nodesCovNorm = nodesCovNorm;
netL.nodesRadAccept_State = nodesRadAccept_State;
netL.nodesRadAccept_Div = nodesRadAccept_Div;
netL.nodesRadAcceptNorm = nodesRadAcceptNorm;
netL.nodesRadAccept = nodesRadAccept;

%%
netL.datanodesNorm = datanodesNorm;
netL.datanodes = datanodes; 
netL.dataColorNode = dataColorNode;
end