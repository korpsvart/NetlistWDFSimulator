function [a] = get_a(elements, b,sampleIndex, inputSignal, a)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
n = size(elements, 1);
% a = zeros(n, 1);
for i=1:n
    type = elements(i, 1);
    switch type
        case 'R'    
            a(i) = 0;
        case 'C'    
            a(i) = b(i);
        case 'L'
            a(i) = -b(i);
        case 'V'
            a(i) = inputSignal(sampleIndex); %small series resistance value
        case 'I'
            a(i) = 10e9*inputSignal(sampleIndex); % large resistance value
    end
end

end

