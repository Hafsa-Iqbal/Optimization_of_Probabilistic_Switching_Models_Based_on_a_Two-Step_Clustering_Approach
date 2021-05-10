function [output] = ReLuF2 (input, variance)

step = 1;
angularCoeff = -1/variance;

if input >= variance
    output = 0;
else
    output = input*angularCoeff + step;
end

end