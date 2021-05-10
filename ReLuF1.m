function [output] = ReLuF1 (input, variance)

step = 0;
angularCoeff = 1/variance;

if input >= variance
    output = 1;
else
    output = input*angularCoeff + step;
end

end