%% GNG  
clc
clear all
close all
% LOAD DATA
load('PMDatafile.mat')

InputData =[structSyncData.Filtered.xPos, structSyncData.Filtered.yPos,...
    structSyncData.Filtered.divxPos, structSyncData.Filtered.divyPos ];

% GNG parameters
params.N = 100;                                                                      %    Number of nodes
params.MaxIt = 10;                                                                   %    Iteration (repetition of input data)
params.L_growing = 1000;                                                             %    Growing rate
params.epsilon_b = 0.05;                                                             %    Movement of winner node
params.epsilon_n = 0.0006;                                                           %    Movement of all other nodes except winner
params.alpha = 0.5;                                                                  %    Decaying global error and utility
params.delta = 0.9995;                                                               %    Decaying local error and utility
params.T = 100;                                                                      %    It could be a function of params.L_growing, e.g., params.LDecay = 2*params.L_growing
params.L_decay = 1000;                                                               %    Decay rate sould be faster than the growing then it will remove extra nodes
params.alpha_utility = 0.0005;                                                       %    It could be a function of params.delta, e.g., params.alpha_utility = 0.9*params.delta
                                                                                     %    Pseudo random shuffling of input data
params.k = 0.7;                                                                      %    Utility threshold  
lossfunction =[];
% GNG processing
for i = 1:1:30
params.seedvector = i; 
netG = GrowingNeuralGasNetwork(InputData, params, true);
F_l = LossFunction(netG.datanodes);
lossfunction = [lossfunction ; F_l ];

[transitionMat,timeMats] = vocabTrans(InputData,netG);
netG.transitionMat = transitionMat;
netG.timeMats = timeMats;
netG.lossfunction = lossfunction;
dataInfo{1,i} = netG; 
i = i+1;

end
Selected_Vocab = dataInfo{1,1};
net = Selected_Vocab;
save('VocabGNG1.mat','net')

save('VocabGNGfull.mat','dataInfo')