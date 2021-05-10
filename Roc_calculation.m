% the Output is ROC, % AUC and ACC
% AUC : Area Under the Curve
% ACC : Accuracy [Roc_f AUC ACC]
function [Roc_f ,AUC, ACC] = Roc_calculation(signal,Ground_truth)
save_TP_f = [];
save_FP_f = [];
condition_positive = sum(Ground_truth);
condition_negative = length(signal) - condition_positive;
for i = 0:0.03:max(signal)
    TP_f = 0;
    FP_f = 0;
    for j = 1: length(signal)
        if signal(j)> i
            if Ground_truth(j)== 1
                TP_f = TP_f + 1;
            else
                FP_f = FP_f +1;
            end
        end
    end
    save_TP_f = [save_TP_f TP_f];
    save_FP_f = [save_FP_f FP_f];
end

true_positive_f  = save_TP_f/condition_positive;
false_positive_f = save_FP_f/condition_negative;
Roc_f = [false_positive_f;true_positive_f];
AUC = -trapz(false_positive_f,true_positive_f);
TN = condition_negative*(1-false_positive_f);
ACC = max(save_TP_f + TN)/(condition_positive + condition_negative);
end

