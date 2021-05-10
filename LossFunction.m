function  F_l = LossFunction(dataNodes) 
        varianceInputData = [];
        N = size(dataNodes,2);                                               % number of nodes at each k
        for i = 1:N                                                          % normalize and ordered data
            tempInputData = dataNodes{1,i};                                  % temprary vector having data samples of a node
            var_InputData = var(tempInputData);                              % variance of first derivative
            varianceInputData =[varianceInputData ; var_InputData];
        end
       
         loss_function = (sigmoid(varianceInputData(:,1)) + (1 - sigmoid(varianceInputData(:,3))) +...
            sigmoid(varianceInputData(:,2)) + (1 - sigmoid(varianceInputData(:,4))) );   % LOSS FUNCTION
        
        avg_loss_function = sum(loss_function)./N;                               % averaging of loss function
        F_l = avg_loss_function;                                                 %  store F_k values in last seed value
       


    end

