clc
clear all
close all
%  UPLOAD DATA
load('PMDatafile.mat')

InputData =[structSyncData.Filtered.xPos, structSyncData.Filtered.yPos,...
    structSyncData.Filtered.divxPos, structSyncData.Filtered.divyPos ];

selectTechnq = 4;  % 1 for GNG , 2 for Kmeans, 3 DBSCAN, 4 SOM



if selectTechnq == 1
    load('VocabGNG1.mat')
    N_opt = [net.nodesMean];
    ColorData = net.dataColorNode;
    NewColorNode = net.dataColorNode;
    dataNodes = net.datanodes;
    Data = net.data;
elseif selectTechnq == 2
    load('VocabKMeans.mat')
    N_opt = netK.nodesMean;
    ColorData = netK.dataColorNode;
    NewColorNode = netK.dataColorNode;
    dataNodes = netK.datanodes;
    Data = netK.InputData;
elseif selectTechnq == 3
    load('VocabDBSCAN.mat')
    N_opt = [net.nodesMean];
    ColorData = net.dataColorNode;
    NewColorNode= net.dataColorNode;
    dataNodes = net.datanodes;
    Data = net.data;
else
    load('VocabSOM.mat')
    N_opt = [net.nodesMean];
    ColorData = net.dataColorNode;
    NewColorNode= net.dataColorNode;
    dataNodes = net.datanodes;
    Data = net.data;
end



%% Hierachical Clustering

HirCluster = linkage(N_opt,'average','chebychev');


figure;
cutoff = median([HirCluster(end-2,3) HirCluster(end-1,3)]);
[H,T,outperm] = dendrogram(HirCluster,'ColorThreshold',cutoff);
hold on
labels_in_level = level_info(HirCluster, size(N_opt,1));

indx1 = 0;

% select clusters at each level
Clusters = cell(1,size(labels_in_level,2)-1);
for i =2:1:size(labels_in_level,2)
    indx2 = size(labels_in_level{i},2);
    indx3 = indx2 + indx1;
    Clusters{1,i-1} = HirCluster(1+indx1:indx3,1:2);
    indx1 = indx1 + indx2;
end

% color data nodes

idx = 0;
for i = 1:1:size(Clusters,2)                                                % to access  each level
    x =  Clusters{1,i};
    for ij = 1:1:size(x,1)                                                  % to acess each cluster who tends to combine
        % data color nodes
        NewColorNode(ColorData == x(ij,1)) = size(dataNodes,2) + idx + ij;
        NewColorNode(ColorData == x(ij,2)) = size(dataNodes,2) + idx + ij;
        ColorData = NewColorNode;
    end
    % data color node
    idx = idx + ij;
    ColorData = NewColorNode;
    ColorDataLevel{1,i} = NewColorNode;
end

% data inside each node
for i =1:1:size(Clusters,2)
    datacolorNode = ColorDataLevel{1,i};
    maxNodes = max(labels_in_level{1, i+1});
    dataNodes = cell(1,maxNodes);
    tempData = Data;
    nData = size(tempData,1);
    for c = 1:nData
        x = datacolorNode(c);
        newSelectedData= tempData(c,:);                                          %   Organize distances between nodes and the first data point in an ascending order
        dataNodes{1,x} = [dataNodes{1,x}; newSelectedData];                      % normalize and ordered data   
    end
    TempDataNodes = dataNodes(~cellfun('isempty',dataNodes));
    dataNodesLevel{1,i} = TempDataNodes;
end


% mean , covariance , radius of each level
for i = 1:1:size(Clusters,2)  % loop on level
    DataNodes = dataNodesLevel{1,i};
    MeanNodes=[];
    CovNodes = cell(1,size(DataNodes,2));
    RadNodes_State =[];
    for xy = 1:1:size(DataNodes,2)
        CurrData = DataNodes{1,xy};
        % mean
        MeanNodes = [MeanNodes; mean(CurrData, 1)];
        % covariance
        CovNodes{1,xy} = cov(CurrData);
        % radius
        RadNodes_State = [RadNodes_State, sqrt(sum((2*std(CurrData(:,1:2))).^2));];
        RadNodes_Div = [RadNodes_State, sqrt(sum((2*std(CurrData(:,1:2))).^2));];
    end
    MeanNodesLevel{1,i} = MeanNodes;
    CovNodesLevel{1,i} = CovNodes;
    RadNodesLevel_State{1,i} = RadNodes_State;
    RadNodesLevel_Div{1,i} = RadNodes_Div;
end

% save hierarchical vocabulary
Hier.Clusters = Clusters;
Hier.DataNodesLevel = dataNodesLevel;
Hier.meanNodesLevel  = MeanNodesLevel;
Hier.covNodesLevel = CovNodesLevel;
Hier.nodesRadLevel_State = RadNodesLevel_State;
Hier.nodesRadLevel_Div = RadNodesLevel_Div;
Hier.ColorDataLevel = ColorDataLevel;
Hier.labels_in_level = labels_in_level;
Hier.HirCluster = HirCluster;
save('VocabHier2.mat','Hier')


% plot results at each level
mycolors = colorcube;
for i = 1:1:size(Clusters,2)
    ColorFig = ColorDataLevel{i};
    MeanFig = MeanNodesLevel{i};
    h = figure;
    hold on
    scatter(Data(:,1),Data(:,2),60,mycolors(ColorFig,:),'.','LineWidth',1)    % colored input data
    scatter(MeanFig(:,1),MeanFig(:,2),250,'+','k','linewidth',2)                     % for the '+' at mean position of nodes
    quiver(MeanFig(:,1),MeanFig(:,2),MeanFig(:,3),MeanFig(:,4),'LineWidth',1.8,'Color','r','AutoScale','on', 'AutoScaleFactor', 0.4)
    grid on
      
end

% variances for each level
InputData = [structSyncData.Filtered.xPos, structSyncData.Filtered.yPos,...
    structSyncData.Filtered.divxPos, structSyncData.Filtered.divyPos ];

%numbers of levels counting the 0
L = size(labels_in_level, 2);
variances = zeros(L, 2);

% f1 and f2 for each level
if selectTechnq == 3 
    InputData = net.nonNoisyData;
end

optvalues = zeros(L, 2);
f = zeros(L, 1);
[optvalues(1, 1),optvalues(1, 2), f(1, 1)] = OptFunctions(net.datanodes, InputData);

for i = 2: L
    [optvalues(i, 1),optvalues(i, 2), f(i, 1)]  = OptFunctions(Hier.DataNodesLevel{1, i-1}, InputData);
end


figure
plot(optvalues(:, 1));
hold on
plot(optvalues(:, 2));
legend 
xlabel('Level');
legend('f1','f2');
set(gca,'XTick',[1 2 3 4 5 6 ] );
set(gca,'XTickLabel',[0 1 2 3 4 5] );

if selectTechnq == 1
    title('Values of f1 and f2 for each level in the GNG case');
elseif selectTechnq == 2
    title('Values of f1 and f2 for each level in the k-means case');
elseif selectTechnq == 3
    title('Values of f1 and f2 for each level in the DBSCAN case');
else
    title('Values of f1 and f2 for each level in the SOM case');
end
